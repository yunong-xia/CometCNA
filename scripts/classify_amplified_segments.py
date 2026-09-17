#!/usr/bin/env python3
"""Reconstruct final CNA segments and classify amplified trajectories."""

import argparse
import csv
import gzip
from collections import defaultdict

CHROMOSOME_LENGTHS = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555,
    5: 181538259, 6: 170805979, 7: 159345973, 8: 145138636,
    9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
    13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345,
    17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}

# open the file in text mode, with gzip support
def open_text(path):
    path = str(path)
    return gzip.open(path, "rt", encoding="ascii") if path.endswith(".gz") else open(path, encoding="ascii")

# read lineages from the lineages table
def read_lineages(path):
    # create a dictionary to store lineages
    lineages = defaultdict(list)
    
    # open the lineages table, given the "path"
    with open_text(path) as handle:
        # each row contains a sample_id, and id. 
        # sample_id is the id of the final sampled cell,
        # id is the id of an ancestor cell in the lineage of the sampled cell.
        for row in csv.DictReader(handle, delimiter="\t"):
            lineages[int(row["sample_id"])].append(int(row["id"]))
    return lineages


# read events from the CNA table
def read_events(path):
    
    # events dictionary
    events = defaultdict(list)
    
    # open CNA table, given the "path"
    with open_text(path) as handle:
        
        # each row contains a CNA event, with the following columns:
        # id: the cell id of the cell that experienced the CNA event
        # cna_event: the type of CNA event (arm_gain, focal_loss, wgd, etc.)
        # chr: the chromosome number of the CNA event
        # arm: the chromosome arm of the CNA event
        # start: the start position of the CNA event
        # end: the end position of the CNA event
        for row in csv.DictReader(handle, delimiter="\t"):
            # event
            event = row["cna_event"]
            # chromosome
            chromosome = int(row["chr"]) if row["chr"] else None
            
            # if event is wgd, append to events with None chromosome and 0 start/end
            if event == "wgd":
                events[int(row["id"])].append((event, None, 0, 0))
            # else
            else:
                # validate chromosome
                if chromosome not in CHROMOSOME_LENGTHS:
                    raise ValueError(f"Unsupported chromosome: {row['chr']}")
                # validate positions
                start = int(row["start"])
                end = int(row["end"])
                if not 0 <= start < end <= CHROMOSOME_LENGTHS[chromosome]:
                    raise ValueError(f"Invalid CNA interval: {row}")
                
                # add event, chromosome, start, end to cell id
                events[int(row["id"])].append((event, chromosome, start, end))
    return events

# classification function
def classify(lineages, events):
    output = []
    for sample_id, child_to_root in sorted(lineages.items()):
        path = list(reversed(child_to_root))
        path_events = [event for cell_id in path for event in events.get(cell_id, [])]
        chromosomes = set(event[1] for event in path_events if event[1] is not None)
        if any(event[0] == "wgd" for event in path_events):
            chromosomes.update(CHROMOSOME_LENGTHS)

        boundaries = {}
        for chromosome in chromosomes:
            boundaries[chromosome] = {0, CHROMOSOME_LENGTHS[chromosome]}
        for _, chromosome, start, end in path_events:
            if chromosome is not None:
                boundaries[chromosome].update((start, end))

        for chromosome in sorted(boundaries):
            points = sorted(boundaries[chromosome])
            for start, end in zip(points, points[1:]):
                trajectory = [2]
                cn = 2
                for event, event_chr, event_start, event_end in path_events:
                    if event == "wgd":
                        cn *= 2
                    elif event_chr == chromosome and event_start <= start and end <= event_end:
                        cn += 1 if event.endswith("gain") else -1
                        cn = max(0, cn)
                    else:
                        continue
                    trajectory.append(cn)
                status = "amplified" if cn > 2 else "non-amplified"
                route_type = "NA"
                if status == "amplified":
                    route_type = "unidirectional" if all(a <= b for a, b in zip(trajectory, trajectory[1:])) else "bidirectional"
                output.append({
                    "sample_id": sample_id, "chr": chromosome, "start": start, "end": end,
                    "cn_trajectory": ",".join(map(str, trajectory)), "final_cn": cn,
                    "amplification_status": status, "route_type": route_type,
                })
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lineages")
    parser.add_argument("cna")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    rows = classify(read_lineages(args.lineages), read_events(args.cna))
    with open(args.output, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t") if rows else None
        if writer:
            writer.writeheader()
            writer.writerows(rows)
    print(f"segments={len(rows)} amplified={sum(r['amplification_status'] == 'amplified' for r in rows)}")


if __name__ == "__main__":
    main()
