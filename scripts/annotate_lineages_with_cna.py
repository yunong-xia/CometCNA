#!/usr/bin/env python3
import argparse
import gzip
import sys


def open_text(path, mode="rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="ascii")
    return open(path, mode, encoding="ascii")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Append CNA event columns to a grouped lineage TSV. The lineage "
            "rows are matched to a CNA TSV/TSV.GZ by cell id."
        )
    )
    parser.add_argument("lineages", help="Grouped lineage TSV/TSV.GZ")
    parser.add_argument("cna", help="CNA TSV/TSV.GZ sorted ascending by id")
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output TSV/TSV.GZ path, or '-' for stdout. Default: stdout",
    )
    return parser.parse_args()


def read_lineages(path):
    with open_text(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_idx = header.index("id")
        except ValueError:
            raise SystemExit("Lineage file is missing required column: id")

        rows = []
        for original_index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= id_idx:
                continue
            rows.append(
                {
                    "original_index": original_index,
                    "id": int(fields[id_idx]),
                    "fields": fields,
                    "events": [],
                }
            )

    return header, rows


class CnaGroupReader:
    def __init__(self, handle, id_idx):
        self.handle = handle
        self.id_idx = id_idx
        self.buffered = None

    def next_group(self):
        first = self.buffered
        self.buffered = None

        while first is None:
            line = self.handle.readline()
            if not line:
                return None, []
            fields = line.rstrip("\n").split("\t")
            if len(fields) > self.id_idx:
                first = fields

        cell_id = int(first[self.id_idx])
        group = [first]

        while True:
            line = self.handle.readline()
            if not line:
                break
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= self.id_idx:
                continue
            next_id = int(fields[self.id_idx])
            if next_id != cell_id:
                self.buffered = fields
                break
            group.append(fields)

        return cell_id, group


def event_summary(cna_rows, cna_header):
    if not cna_rows:
        return ["0", "", "", "", "", ""]

    event_idx = cna_header.index("cna_event")
    chr_idx = cna_header.index("chr")
    arm_idx = cna_header.index("arm")
    start_idx = cna_header.index("start")
    end_idx = cna_header.index("end")

    return [
        str(len(cna_rows)),
        ";".join(row[event_idx] for row in cna_rows),
        ";".join(row[chr_idx] for row in cna_rows),
        ";".join(row[arm_idx] for row in cna_rows),
        ";".join(row[start_idx] for row in cna_rows),
        ";".join(row[end_idx] for row in cna_rows),
    ]


def annotate_with_cna(rows, cna_path):
    rows_by_id = sorted(rows, key=lambda row: row["id"])
    row_pos = 0
    matched_rows = 0
    matched_events = 0
    cna_header = []

    with open_text(cna_path, "rt") as handle:
        cna_header = handle.readline().rstrip("\n").split("\t")
        required = ["id", "cna_event", "chr", "arm", "start", "end"]
        missing = [name for name in required if name not in cna_header]
        if missing:
            raise SystemExit(f"CNA file is missing required columns: {','.join(missing)}")
        cna_id_idx = cna_header.index("id")

        cna_reader = CnaGroupReader(handle, cna_id_idx)
        cna_id, cna_group = cna_reader.next_group()

        while row_pos < len(rows_by_id):
            lineage_id = rows_by_id[row_pos]["id"]

            if cna_id is None:
                row_pos += 1
                continue

            if lineage_id < cna_id:
                row_pos += 1
                continue

            if lineage_id > cna_id:
                cna_id, cna_group = cna_reader.next_group()
                continue

            window_start = row_pos
            while row_pos < len(rows_by_id) and rows_by_id[row_pos]["id"] == lineage_id:
                row_pos += 1

            for row in rows_by_id[window_start:row_pos]:
                row["events"] = cna_group
                matched_rows += 1
            matched_events += len(cna_group)

            cna_id, cna_group = cna_reader.next_group()

    return matched_rows, matched_events, cna_header


def write_output(path, header, rows, cna_header):
    extra_header = [
        "cna_event_count",
        "cna_events",
        "cna_chrs",
        "cna_arms",
        "cna_starts",
        "cna_ends",
    ]

    with open_text(path, "wt") as handle:
        handle.write("\t".join(header + extra_header) + "\n")
        for row in sorted(rows, key=lambda item: item["original_index"]):
            handle.write(
                "\t".join(row["fields"] + event_summary(row["events"], cna_header))
                + "\n"
            )


def main():
    args = parse_args()
    header, rows = read_lineages(args.lineages)
    matched_rows, matched_events, cna_header = annotate_with_cna(rows, args.cna)
    write_output(args.output, header, rows, cna_header)

    print(
        f"lineage_rows={len(rows)} rows_with_cna={matched_rows} "
        f"cna_events_matched={matched_events}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
