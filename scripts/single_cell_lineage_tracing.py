#!/usr/bin/env python3
"""single cell lineage tracing."""
import argparse
import csv
import gzip
import pandas as pd
from collections import defaultdict
from pathlib import Path
from bisect import bisect_left
import numpy as np

DEFAULT_CYTOBAND_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "cytoband_hg38.txt"
)
# combine population and dead cells data into a single dictionary

def sample_extant_cells(population_path, sample_size, rng):
    sampled_cells = []
    n_extant_cells_encountered = 0
    opener = gzip.open if str(population_path).endswith('.gz') else open
    with opener(population_path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            if float(row['death']) != 0:
                continue  # skip if a dead cell is encountered

            n_extant_cells_encountered += 1

            # reservoir sampling algorithm to sample cells uniformly at random
            # put the first sample_size cells into the sampled_cells list
            if len(sampled_cells) < sample_size:
                sampled_cells.append(int(row['id']))
            # afterwards, use the val as a random index to replace a cell in the sampled_cells list with probability sample_size/n_extant_cells_encountered
            else:
                val = rng.integers(1, n_extant_cells_encountered + 1)
                # only replace a cell in the sampled_cells list if val <= sample_size
                if val <= sample_size:
                    # Replace a random cell in the sampled_cells list
                    sampled_cells[val - 1] = int(row['id'])

    return sampled_cells

def load_population_and_dead_cells(population_path, dead_cells_path):
    # Load population data
    opener = gzip.open if str(population_path).endswith('.gz') else open
    with opener(population_path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        population_data = {int(row['id']): row for row in reader}

    # Load dead cells data
    with opener(dead_cells_path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        dead_cells_data = {int(row['id']): row for row in reader}

    all_cells = population_data.copy()
    all_cells.update(dead_cells_data)


    # dictionary data structure already remove duplicates between population and dead cells
    cell_parent_dict = {
        cell_id: int(row['ancestor'])
        for cell_id, row in all_cells.items()
    }

    return cell_parent_dict

def load_cna_events(cna_events_path):
    # Load CNA events data
    opener = gzip.open if str(cna_events_path).endswith('.gz') else open
    with opener(cna_events_path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        cna_events_data = {int(row['id']): row for row in reader}
    return cna_events_data

# from ancestor to descendant, trace the lineage of a given cell
def trace_lineage(cell_parent_dict, cell_id):
    lineage = []
    current_cell = cell_id
    while current_cell in cell_parent_dict:
        lineage.append(current_cell)
        current_cell = cell_parent_dict[current_cell]
    lineage.append(current_cell)  # Add the root ancestor
    return lineage[::-1]  # Reverse to get from ancestor to descendant

#  CNA events along the lineage from ancestor to descendant
def trace_cna_events(lineage, cna_events_dict):
    lineage_cna_events = []
    for cell_id in lineage:
        if cell_id in cna_events_dict:
            lineage_cna_events.append(cna_events_dict[cell_id])
    return lineage_cna_events


def initialize_diploid_genome(cytoband_path):
    arm_bounds = {}

    with open(cytoband_path, newline='') as handle:
        for chrom, start, end, band in csv.reader(handle, delimiter='\t'):
            if chrom not in {f'chr{i}' for i in range(1, 23)}:
                continue

            arm = chrom[3:] + band[0]
            start, end = int(start), int(end)

            if arm not in arm_bounds:
                arm_bounds[arm] = [start, end]
            else:
                arm_bounds[arm][0] = min(arm_bounds[arm][0], start)
                arm_bounds[arm][1] = max(arm_bounds[arm][1], end)

    bps = {arm: [start, end] for arm, (start, end) in arm_bounds.items()}
    cn_states = {arm: [2] for arm in arm_bounds}

    return cn_states, bps

def split_at(arm, position, cn_states, bps):
    boundaries = bps[arm]
    index = bisect_left(boundaries, position)

    if index < len(boundaries) and boundaries[index] == position:
        return

    if 0 < index < len(boundaries):
        boundaries.insert(index, position)
        cn_states[arm].insert(index, cn_states[arm][index - 1])


def apply_delta(arm, start, end, delta, cn_states, bps):
    split_at(arm, start, cn_states, bps)
    split_at(arm, end, cn_states, bps)

    boundaries = bps[arm]
    states = cn_states[arm]

    for i, cn in enumerate(states):
        if boundaries[i] >= start and boundaries[i + 1] <= end and cn > 0:
            states[i] = max(0, cn + delta)


def reconstruct_copy_numbers(lineage_cna_events, cytoband_path=DEFAULT_CYTOBAND_PATH):
    cn_states, bps = initialize_diploid_genome(cytoband_path)

    for event in lineage_cna_events:
        event_type = event['cna_event']

        if event_type == 'wgd':
            for arm in cn_states:
                cn_states[arm] = [cn * 2 for cn in cn_states[arm]]
            continue

        if event_type.endswith('_gain'):
            delta = 1
        elif event_type.endswith('_loss'):
            delta = -1
        else:
            raise ValueError(f'Unknown CNA event: {event_type}')

        if event_type.startswith('chromosome_'):
            chromosome = str(int(event['chr']))
            affected_arms = [chromosome + 'p', chromosome + 'q']
        else:
            affected_arms = [event['arm']]

        start = int(event['start'])
        end = int(event['end'])

        for arm in affected_arms:
            apply_delta(arm, start, end, delta, cn_states, bps)

    return cn_states, bps

def copy_number_table(cell_id, cn_states, bps):
    rows = []

    for arm in sorted(cn_states, key=lambda name: (int(name[:-1]), name[-1])):
        boundaries = bps[arm]
        states = cn_states[arm]

        if len(boundaries) != len(states) + 1:
            raise ValueError(f'Invalid breakpoints/CN states for {arm}')

        for start, end, cn in zip(boundaries[:-1], boundaries[1:], states):
            rows.append({
                'cell_id': cell_id,
                'chr': int(arm[:-1]),
                'arm': arm,
                'start': start,
                'end': end,
                'cn': cn,
            })

    return pd.DataFrame(rows, columns=['cell_id', 'chr', 'arm', 'start', 'end', 'cn'])

def main():
    # runs/tumor_14q_selection_coef_0.1_N1e5_seed100001

    argument_parser = argparse.ArgumentParser(description="Single cell lineage tracing.")
    argument_parser.add_argument("--population_path", type=str, required=True, help="Path to the population.tsv.gz file.")
    argument_parser.add_argument("--dead_cells_path", type=str, required=True, help="Path to the dead_cells.tsv.gz file.")
    argument_parser.add_argument("--cna_events_path", type=str, required=True, help="Path to the cna.tsv.gz file.")
    argument_parser.add_argument("--sample_size", type=int, default=10, help="Number of extant cells to sample.")
    argument_parser.add_argument("--output_path", type=str, required=True, help="Path to save the sampled cells copy number table.")
    args = argument_parser.parse_args()
    if args.sample_size < 1:
        argument_parser.error("--sample_size must be at least 1")

    population_path = Path(args.population_path)
    dead_cells_path = Path(args.dead_cells_path)

     # sample extant cells
    sampled_cells = sample_extant_cells(population_path, args.sample_size, np.random.default_rng(42))
    if len(sampled_cells) != args.sample_size:
        argument_parser.error(
            f"Requested {args.sample_size} cells, but found only {len(sampled_cells)} extant cells"
        )
    print(f"Sampled cells: {sampled_cells}")

    cell_parent_dict = load_population_and_dead_cells(population_path, dead_cells_path)

    cna_events_path = Path(args.cna_events_path)
    cna_events_dict = load_cna_events(cna_events_path)




    table = pd.DataFrame(columns=['cell_id', 'chr', 'arm', 'start', 'end', 'cn'])
    for cell_id in sampled_cells:
        print(f"Processing cell {cell_id}...")
        cn_states, bps = reconstruct_copy_numbers(trace_cna_events(trace_lineage(cell_parent_dict, cell_id), cna_events_dict))
        table = pd.concat([table, copy_number_table(cell_id, cn_states, bps)], ignore_index=True)
    table.to_csv(args.output_path, sep="\t", index=False)
    return


if __name__ == "__main__":
    main()
