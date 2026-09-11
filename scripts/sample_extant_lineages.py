#!/usr/bin/env python3
import argparse
import csv
import gzip
import random
import sys
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii")
    return open(path, "rt", encoding="ascii")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Sample extant cells from a population file, then trace each "
            "sampled cell's lineage through a cell table sorted by id descending."
        )
    )
    parser.add_argument("population", help="Population TSV/TSV.GZ file")
    parser.add_argument(
        "sorted_cells",
        help="TSV/TSV.GZ with id in column 4, sorted descending by id",
    )
    parser.add_argument("sample_size", type=int, nargs="?", help="Number of extant cells to sample")
    parser.add_argument("--sample-ids", help="TSV containing an id column, instead of random sampling")
    parser.add_argument("--cell-id", type=int, help="Trace one cell ID, including a dead ancestor")
    parser.add_argument("--mrca-output", help="Write the shared most recent ancestor as a TSV")
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output TSV path, or '-' for stdout. Default: stdout",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible sampling",
    )
    return parser.parse_args()


def reservoir_sample_extant(population_path, sample_size, rng):
    reservoir = []
    extant_seen = 0

    with open_text(population_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_idx = header.index("id")
            death_idx = header.index("death")
        except ValueError as exc:
            raise SystemExit(f"Missing required column in population header: {exc}")

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(id_idx, death_idx):
                continue
            if fields[death_idx] != "0":
                continue

            extant_seen += 1
            cell_id = fields[id_idx]
            if len(reservoir) < sample_size:
                reservoir.append(cell_id)
            else:
                j = rng.randrange(extant_seen)
                if j < sample_size:
                    reservoir[j] = cell_id

    if extant_seen < sample_size:
        raise SystemExit(
            f"Requested {sample_size} extant cells, but only found {extant_seen}"
        )

    return reservoir, extant_seen


def read_sample_ids(sample_path, population_path):
    with open_text(sample_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "id" not in (reader.fieldnames or []):
            raise SystemExit("Sample table must contain an id column")
        ids = [str(int(row["id"])) for row in reader]
    if not ids or len(set(ids)) != len(ids) or any(int(x) <= 0 for x in ids):
        raise SystemExit("Sample IDs must be nonempty, positive, and unique")
    missing = set(ids)
    with open_text(population_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"id", "death"}.issubset(reader.fieldnames or []):
            raise SystemExit("Population must contain id and death columns")
        for row in reader:
            if row["id"] in missing and float(row["death"]) == 0:
                missing.remove(row["id"])
    if missing:
        raise SystemExit(f"Sample IDs not found as extant cells: {sorted(missing, key=int)[:10]}")
    return ids


def read_cell_id(cell_id, population_path):
    if cell_id <= 0:
        raise SystemExit("cell ID must be positive")
    with open_text(population_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "id" not in (reader.fieldnames or []):
            raise SystemExit("Population must contain an id column")
        for row in reader:
            if row["id"] == str(cell_id):
                return [str(cell_id)]
    raise SystemExit(f"Cell ID not found in population: {cell_id}")


def trace_lineages(sorted_cells_path, sample_ids, output_handle, mrca_handle=None):
    targets = defaultdict(set)
    for sample_id in sample_ids:
        targets[sample_id].add(sample_id)

    rows_written = 0
    lineages = {sample_id: [] for sample_id in sample_ids}

    with open_text(sorted_cells_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_idx = header.index("id")
            ancestor_idx = header.index("ancestor")
        except ValueError as exc:
            raise SystemExit(f"Missing required column in sorted-cells header: {exc}")

        previous_id = float("inf")
        for line in handle:
            if not targets:
                break

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(id_idx, ancestor_idx):
                continue

            cell_id = fields[id_idx]
            if int(cell_id) > previous_id:
                raise SystemExit("Cell table must be sorted by id descending")
            previous_id = int(cell_id)
            if cell_id not in targets:
                continue

            sample_set = targets.pop(cell_id)
            ancestor = fields[ancestor_idx]
            if not 0 <= int(ancestor) < int(cell_id):
                raise SystemExit(f"Invalid parent ID for cell {cell_id}: {ancestor}")

            for sample_id in sorted(sample_set, key=int):
                lineages[sample_id].append(fields)
                rows_written += 1

            if ancestor != "0":
                targets[ancestor].update(sample_set)

        if mrca_handle is not None:
            if targets:
                raise SystemExit("Incomplete lineages: cannot determine MRCA")
            common = None
            for lineage in lineages.values():
                ids = {fields[id_idx] for fields in lineage}
                common = ids if common is None else common & ids
            if not common:
                raise SystemExit("No shared cell ancestor (ID 0 is a sentinel, not a cell)")
            # IDs increase along descent, so the largest shared ID is the MRCA.
            mrca_id = max(common, key=int)
            mrca_row = next(row for row in lineages[sample_ids[0]] if row[id_idx] == mrca_id)
            mrca_handle.write("sample_count\t" + "\t".join(header) + "\n")
            mrca_handle.write(str(len(sample_ids)) + "\t" + "\t".join(mrca_row) + "\n")

        output_handle.write("sample_id\tlineage_index\t" + "\t".join(header) + "\n")
        for sample_id in sample_ids:
            lineage = lineages[sample_id]
            lineage.sort(key=lambda fields: int(fields[id_idx]), reverse=True)
            for lineage_index, fields in enumerate(lineage):
                output_handle.write(
                    sample_id
                    + "\t"
                    + str(lineage_index)
                    + "\t"
                    + "\t".join(fields)
                    + "\n"
                )

    return rows_written, targets


def main():
    args = parse_args()
    selected_modes = sum(bool(value) for value in (args.sample_ids, args.cell_id is not None))
    if selected_modes > 1 or (selected_modes and args.sample_size is not None):
        raise SystemExit("Use only one of sample_size, --sample-ids, or --cell-id")
    if not selected_modes and (args.sample_size is None or args.sample_size <= 0):
        raise SystemExit("sample_size must be positive")

    rng = random.Random(args.seed)
    if args.cell_id is not None:
        sample_ids = read_cell_id(args.cell_id, args.population)
        extant_seen = "not_counted"
    elif args.sample_ids:
        sample_ids = read_sample_ids(args.sample_ids, args.population)
        extant_seen = "not_counted"
    else:
        sample_ids, extant_seen = reservoir_sample_extant(
            args.population, args.sample_size, rng
        )

    from contextlib import ExitStack
    with ExitStack() as stack:
        output_handle = sys.stdout if args.output == "-" else stack.enter_context(open(args.output, "w", encoding="ascii"))
        mrca_handle = stack.enter_context(open(args.mrca_output, "w", encoding="ascii")) if args.mrca_output else None
        rows_written, unresolved = trace_lineages(
            args.sorted_cells, sample_ids, output_handle, mrca_handle
        )

    print(
        f"sampled_extant={len(sample_ids)} total_extant_seen={extant_seen} "
        f"lineage_rows_written={rows_written} unresolved_targets={len(unresolved)}",
        file=sys.stderr,
    )
    if unresolved:
        preview = ",".join(list(unresolved.keys())[:10])
        print(f"unresolved target ids: {preview}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
