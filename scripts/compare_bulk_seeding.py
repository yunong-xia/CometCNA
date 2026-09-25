#!/usr/bin/env python3
"""Compare a bulk sample's MRCA with seeding cells using recorded ancestry."""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def read_table(path, required):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: required columns: {required}')
        yield from reader


def read_ids(path):
    ids = list(dict.fromkeys(int(row['id']) for row in read_table(path, ['id'])))
    if not ids or min(ids) <= 0:
        raise ValueError(f'{path}: expected positive cell IDs')
    return ids


def load_parents(paths):
    parents = {}
    for path in paths:
        for row in read_table(path, ['id', 'ancestor']):
            cell, parent = int(row['id']), int(row['ancestor'])
            if not 0 <= parent < cell:
                raise ValueError(f'Invalid ancestry: {cell} -> {parent}')
            if cell in parents and parents[cell] != parent:
                raise ValueError(f'Conflicting parents for cell {cell}')
            parents[cell] = parent
    return parents


def lineage(cell, parents):
    result = []
    while cell:
        result.append(cell)
        if cell not in parents:
            raise ValueError(f'Missing ancestry for cell {cell}; include dead_cells table')
        cell = parents[cell]
    return result


def common_ancestor(paths):
    common = set(paths[0])
    for path in paths[1:]:
        common.intersection_update(path)
    if not common:
        raise ValueError('No shared cell ancestor; ID 0 is a sentinel')
    return next(cell for cell in paths[0] if cell in common)


def write_table(path, header, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


def run(args):
    bulk_ids, seed_ids = read_ids(args.bulk_ids), read_ids(args.seeding_ids)
    parents = load_parents([args.population, args.dead_cells])
    missing = set(bulk_ids)
    for row in read_table(args.population, ['id', 'death']):
        if float(row['death']) == 0:
            missing.discard(int(row['id']))
    if missing:
        raise ValueError(f'Bulk IDs are not extant in population: {sorted(missing)[:10]}')
    # Intersect one path at a time; do not retain all bulk-cell lineages.
    first = lineage(bulk_ids[0], parents)
    common = set(first)
    for cell in bulk_ids[1:]:
        common.intersection_update(lineage(cell, parents))
    if not common:
        raise ValueError('Bulk sample has no shared cell ancestor')
    bulk_mrca = next(cell for cell in first if cell in common)
    bulk_path = lineage(bulk_mrca, parents)
    seed_paths = {cell: lineage(cell, parents) for cell in seed_ids}
    pairs, branches = [], []
    for seed, seed_path in seed_paths.items():
        shared = common_ancestor([bulk_path, seed_path])
        primary = list(reversed(bulk_path[:bulk_path.index(shared)]))
        seeding = list(reversed(seed_path[:seed_path.index(shared)]))
        pairs.append((seed, bulk_mrca, shared, len(primary), len(seeding)))
        for branch, path in [('primary', primary), ('seeding', seeding)]:
            for step, cell in enumerate(path, 1):
                branches.append((seed, bulk_mrca, shared, branch, step, cell, parents[cell]))
    relevant = {row[5] for row in branches}
    events = defaultdict(list)
    cna_columns = ['id', 'cna_event', 'chr', 'arm', 'start', 'end']
    for row in read_table(args.cna, cna_columns):
        cell = int(row['id'])
        if cell in relevant:
            events[cell].append([row[col] for col in cna_columns[1:]])
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / 'bulk_mrca.tsv', ['sample_count', 'id'], [(len(bulk_ids), bulk_mrca)])
    lineage_header = ['sample_id', 'lineage_index', 'id', 'ancestor']
    write_table(out / 'bulk_mrca_lineage.tsv', lineage_header,
                ((bulk_mrca, i, cell, parents[cell]) for i, cell in enumerate(bulk_path)))
    write_table(out / 'seeding_lineages.tsv', lineage_header,
                ((seed, i, cell, parents[cell]) for seed, path in seed_paths.items()
                 for i, cell in enumerate(path)))
    write_table(out / 'pair_mrcas.tsv',
                ['seeding_id', 'bulk_mrca_id', 'primary_metastasis_mrca_id',
                 'primary_branch_cells', 'seeding_branch_cells'], pairs)
    branch_header = ['seeding_id', 'bulk_mrca_id', 'primary_metastasis_mrca_id',
                     'branch', 'branch_step', 'id', 'ancestor']
    write_table(out / 'branches.tsv', branch_header, branches)
    write_table(out / 'branch_cna_events.tsv', branch_header + cna_columns[1:],
                (list(row) + event for row in branches for event in events[row[5]]))
    print(f'bulk_mrca={bulk_mrca} seeding_cells={len(seed_ids)} output={out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['population', 'dead-cells', 'bulk-ids', 'seeding-ids', 'cna', 'outdir']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
