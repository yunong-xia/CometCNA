#!/usr/bin/env python3
"""single cell lineage tracing."""
import argparse
import csv
import gzip
import pandas as pd
from collections import defaultdict
from pathlib import Path
from bisect import bisect_left, bisect_right
import numpy as np
from copy import deepcopy

import single_cell_lineage_tracing as sclt

DEFAULT_CYTOBAND_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "cytoband_hg38.txt"
)


def load_sampled_seeding_cells(sampled_seeding_cells_path):
    opener = gzip.open if str(sampled_seeding_cells_path).endswith('.gz') else open
    seeding_cells_id = []
    with opener(sampled_seeding_cells_path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        seeding_cells_id = [int(row['id']) for row in reader]
    return seeding_cells_id

def trace_subtrajectory_lineages(cell_parent_dict, primary_cell, seeding_cell):
    primary_cell_lineage = sclt.trace_lineage(cell_parent_dict, primary_cell)
    seeding_cell_lineage = sclt.trace_lineage(cell_parent_dict, seeding_cell)

    # truncal lineage
    truncal_lineage = sorted(set(primary_cell_lineage) & set(seeding_cell_lineage))

    # subclonal lineage to the primary cell and the seeding cell
    primary_subclonal_lineage = sorted(set(primary_cell_lineage) - set(seeding_cell_lineage))
    seeding_subclonal_lineage = sorted(set(seeding_cell_lineage) - set(primary_cell_lineage))

    return truncal_lineage, primary_subclonal_lineage, seeding_subclonal_lineage

def trace_subtrajectory_cna_events(cell_parent_dict, cna_events_dict, primary_cell, seeding_cell):
    truncal_lineage, primary_subclonal_lineage, seeding_subclonal_lineage = trace_subtrajectory_lineages(cell_parent_dict, primary_cell, seeding_cell)
    # truncal CNA events
    truncal_cna_events = sclt.trace_cna_events(truncal_lineage, cna_events_dict)

    # trace the CNA events along each lineage
    primary_subclonal_cna_events = sclt.trace_cna_events(primary_subclonal_lineage, cna_events_dict)
    seeding_subclonal_cna_events = sclt.trace_cna_events(seeding_subclonal_lineage, cna_events_dict)

    return truncal_cna_events, primary_subclonal_cna_events, seeding_subclonal_cna_events


# trace the lineage, and get the cn states, bps, and trajectories for primary vs seeding
def get_primary_vs_seeding_genome_wide_cn_trajectories(cell_parent_dict, cna_events_dict, primary_cell, seeding_cell):
    truncal_cna_events, primary_subclonal_cna_events, seeding_subclonal_cna_events = trace_subtrajectory_cna_events(cell_parent_dict, cna_events_dict, primary_cell, seeding_cell)
    truncal_cn_states, truncal_bps, truncal_trajectory_types = sclt.reconstruct_copy_numbers(truncal_cna_events)
    
    
    
    primary_cn_states, primary_bps, primary_subclonal_trajectory_types = sclt.reconstruct_copy_numbers(primary_subclonal_cna_events,
                                                                                                       starting_cn_states=deepcopy(truncal_cn_states),
                                                                                                       starting_bps=deepcopy(truncal_bps))



    seeding_cn_states, seeding_bps, seeding_subclonal_trajectory_types = sclt.reconstruct_copy_numbers(seeding_subclonal_cna_events,
                                                                                                    starting_cn_states=deepcopy(truncal_cn_states),
                                                                                                    starting_bps=deepcopy(truncal_bps))


    combined_bps = {}
    combined_types = {}
    primary_cn_states_on_union_bps = {}
    seeding_cn_states_on_union_bps = {}
    truncal_cn_states_on_union_bps = {}
    
    primary_branch_types = {}
    seeding_branch_types = {}
    truncal_types = {}

    for arm in primary_bps:

        # take unions of breakpoints
        boundaries = sorted(set(primary_bps[arm]) | set(seeding_bps[arm]) | set(truncal_bps[arm]))

        # initialization
        labels = []
        arm_primary_branch_types = []
        arm_seeding_branch_types = []
        arm_truncal_types = []
        
        arm_seeding_cn_states_on_union_bps = []
        arm_truncal_cn_states_on_union_bps = []
        
        arm_primary_cn_states_on_union_bps = []
        arm_seeding_cn_states_on_union_bps = []
        arm_truncal_cn_states_on_union_bps = []

        # for each pair of consecutive union breakpoints
        for start, end in zip(boundaries[:-1], boundaries[1:]):

            # index of primary cell bp <= the union bp (indexing bps)
            p_i = bisect_right(primary_bps[arm], start) - 1
            # index of seeding cell bp <= the union bp (indexing bps)
            s_i = bisect_right(seeding_bps[arm], start) - 1
            # index of the MRCA bp <= the union bp (indexing bps)
            t_i = bisect_right(truncal_bps[arm], start) - 1

            is_bdt = (
                primary_subclonal_trajectory_types[arm][p_i] == "BDT"
                or seeding_subclonal_trajectory_types[arm][s_i] == "BDT"
            )
            labels.append("BDT" if is_bdt else "UDT")
            
            arm_primary_branch_types.append(primary_subclonal_trajectory_types[arm][p_i])
            arm_seeding_branch_types.append(seeding_subclonal_trajectory_types[arm][s_i])
            arm_truncal_types.append(truncal_trajectory_types[arm][t_i])
            
            arm_primary_cn_states_on_union_bps.append(primary_cn_states[arm][p_i])
            arm_seeding_cn_states_on_union_bps.append(seeding_cn_states[arm][s_i])
            arm_truncal_cn_states_on_union_bps.append(truncal_cn_states[arm][t_i])
            
        combined_bps[arm] = boundaries
        combined_types[arm] = labels
        
        truncal_types[arm] = arm_truncal_types
        primary_subclonal_trajectory_types[arm] = arm_primary_branch_types
        seeding_subclonal_trajectory_types[arm] = arm_seeding_branch_types
        
        primary_cn_states_on_union_bps[arm] = arm_primary_cn_states_on_union_bps
        seeding_cn_states_on_union_bps[arm] = arm_seeding_cn_states_on_union_bps
        truncal_cn_states_on_union_bps[arm] = arm_truncal_cn_states_on_union_bps
        

    primary_vs_seeding_results = {
        "primary_cell": primary_cell,
        "seeding_cell": seeding_cell,
        "primary_cn_states_on_union_bps": primary_cn_states_on_union_bps,
        "seeding_cn_states_on_union_bps": seeding_cn_states_on_union_bps,
        "truncal_cn_states_on_union_bps": truncal_cn_states_on_union_bps,
        "primary_subclonal_trajectory_types": primary_subclonal_trajectory_types,
        "seeding_subclonal_trajectory_types": seeding_subclonal_trajectory_types,
        "truncal_types": truncal_types,
        "combined_bps": combined_bps,
        "combined_types": combined_types
    }
    
    return primary_vs_seeding_results

def primary_vs_seeding_table(primary_vs_seeding_results):
    
    primary_cell = primary_vs_seeding_results['primary_cell']
    seeding_cell = primary_vs_seeding_results['seeding_cell']
    
    combined_bps = primary_vs_seeding_results['combined_bps']
    combined_types = primary_vs_seeding_results['combined_types']
    
    primary_cn_states_on_union_bps = primary_vs_seeding_results['primary_cn_states_on_union_bps']
    seeding_cn_states_on_union_bps = primary_vs_seeding_results['seeding_cn_states_on_union_bps']
    truncal_cn_states_on_union_bps = primary_vs_seeding_results['truncal_cn_states_on_union_bps']
    
    primary_subclonal_trajectory_types = primary_vs_seeding_results['primary_subclonal_trajectory_types']
    seeding_subclonal_trajectory_types = primary_vs_seeding_results['seeding_subclonal_trajectory_types']
    truncal_types = primary_vs_seeding_results['truncal_types']
    
    
    rows = []
    
    for arm in sorted(combined_types, key=lambda name: (int(name[:-1]), name[-1])):
        boundaries = combined_bps[arm]
        types = combined_types[arm]
        
        arm_primary_cn_states_on_union_bps = primary_cn_states_on_union_bps[arm]
        arm_seeding_cn_states_on_union_bps = seeding_cn_states_on_union_bps[arm]
        arm_truncal_cn_states_on_union_bps = truncal_cn_states_on_union_bps[arm]
        
        arm_primary_subclonal_trajectory_types = primary_subclonal_trajectory_types[arm]
        arm_seeding_subclonal_trajectory_types = seeding_subclonal_trajectory_types[arm]
        arm_truncal_types = truncal_types[arm]

        if len(boundaries) != len(types) + 1:
            raise ValueError(f'Invalid breakpoints/CN states for {arm}')

        for start, end, type, primary_cn, seeding_cn, truncal_cn, primary_type, seeding_type, truncal_type in zip(boundaries[:-1], 
                                                            boundaries[1:], 
                                                            types, 
                                                            arm_primary_cn_states_on_union_bps, 
                                                            arm_seeding_cn_states_on_union_bps,
                                                            arm_truncal_cn_states_on_union_bps,
                                                            arm_primary_subclonal_trajectory_types,
                                                            arm_seeding_subclonal_trajectory_types,
                                                            arm_truncal_types):
            rows.append({
                'primary_cell': primary_cell,
                'seeding_cell': seeding_cell,
                'chr': int(arm[:-1]),
                'arm': arm,
                'start': start,
                'end': end,
                'primary_cell_cn': primary_cn,
                'seeding_cell_cn': seeding_cn,
                'truncal_cell_cn': truncal_cn,
                'primary_subtrajectory_type': primary_type,
                'seeding_subtrajectory_type': seeding_type,
                'truncal_subtrajectory_type': truncal_type,
                'primary_vs_seeding_subclonal_trajectory_type': type
            })

    return pd.DataFrame(rows, columns=['primary_cell', 
                                       'seeding_cell', 
                                       'chr', 
                                       'arm', 
                                       'start', 
                                       'end', 
                                       'primary_cell_cn', 
                                       'seeding_cell_cn', 
                                       'truncal_cell_cn',
                                       'primary_subtrajectory_type',
                                       'seeding_subtrajectory_type',
                                       'truncal_subtrajectory_type',
                                       'primary_vs_seeding_subclonal_trajectory_type'])

def main():
    argument_parser = argparse.ArgumentParser(description="Single cell lineage tracing.")
    argument_parser.add_argument("--population_path", type=str, required=True, help="Path to the population.tsv.gz file.")
    argument_parser.add_argument("--dead_cells_path", type=str, required=True, help="Path to the dead_cells.tsv.gz file.")
    argument_parser.add_argument("--sampled_seeding_cells_path", type=str, required=True, help="Path to the sampled seeding cells file.")
    argument_parser.add_argument("--cna_events_path", type=str, required=True, help="Path to the cna.tsv.gz file.")
    argument_parser.add_argument("--sample_size", type=int, default=10, help="Number of extant cells to sample.")
    argument_parser.add_argument("--output_path", type=str, required=True, help="Path to save the sampled cells copy number table.")
    argument_parser.add_argument("--seed", type=int, default=None, help="Random seed for sampling.")
    args = argument_parser.parse_args()
    if args.sample_size < 1:
        argument_parser.error("--sample_size must be at least 1")


    if args.seed is not None:
        np.random.seed(args.seed)
    else:
        np.random.seed()


    population_path = Path(args.population_path)
    dead_cells_path = Path(args.dead_cells_path)

    # sample extant cells
    primary_sampled_cells = sclt.sample_extant_cells(population_path, args.sample_size, np.random.default_rng(args.seed))
    if len(primary_sampled_cells) != args.sample_size:
        argument_parser.error(
            f"Requested {args.sample_size} cells, but found only {len(primary_sampled_cells)} extant cells"
        )
    print(f"Sampled cells: {primary_sampled_cells}")

    # sampled seeding cells.
    sampled_seeding_cells_id = load_sampled_seeding_cells(args.sampled_seeding_cells_path)
    sampled_seeding_cells = sampled_seeding_cells_id

    # (primary, seeding) pairs
    primary_seeding_pairs = [(primary, seeding) for primary in primary_sampled_cells for seeding in sampled_seeding_cells]
    print(f"Number of primary-Seeding pairs: {len(primary_seeding_pairs)}")


    # read parent-child relationships
    cell_parent_dict = sclt.load_population_and_dead_cells(population_path, dead_cells_path)

    # read CNA events
    cna_events_path = Path(args.cna_events_path)
    cna_events_dict = sclt.load_cna_events(cna_events_path)

    # trace subclonal CNA events
    tables = []
    for pair in primary_seeding_pairs:
        primary_cell, seeding_cell = pair
        primary_vs_seeding_results = get_primary_vs_seeding_genome_wide_cn_trajectories(cell_parent_dict, cna_events_dict,
                                               primary_cell, seeding_cell)
        tables.append(primary_vs_seeding_table(primary_vs_seeding_results))
    columns = [
        'primary_cell', 
        'seeding_cell', 
        'chr', 
        'arm', 
        'start', 
        'end', 
        'primary_cell_cn', 
        'seeding_cell_cn', 
        'truncal_cell_cn',
        'primary_subtrajectory_type',
        'seeding_subtrajectory_type',
        'truncal_subtrajectory_type',
        'primary_vs_seeding_subclonal_trajectory_type'
    ]
    print(tables)
    table = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=columns)
    table.to_csv(args.output_path, sep="\t", index=False)

    return

if __name__ == "__main__":
    main()
