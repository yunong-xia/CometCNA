#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: sort_cells_by_id_desc.sh INPUT.tsv.gz [OUTPUT.tsv.gz]

Sort a simulator TSV/TSV.GZ file by cell id in descending numeric order.
The input header is preserved. The cell id column is expected to be column 4.

Environment variables:
  SORT_MEM  Memory budget for sort, default: 4G
  TMPDIR    Temporary directory for sort spill files, default: /tmp
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

input=$1
output=${2:-}

if [[ ! -f "$input" ]]; then
  echo "Input file not found: $input" >&2
  exit 1
fi

if [[ -z "$output" ]]; then
  if [[ "$input" == *.tsv.gz ]]; then
    output="${input%.tsv.gz}.sorted.desc.tsv.gz"
  elif [[ "$input" == *.gz ]]; then
    output="${input%.gz}.sorted.desc.tsv.gz"
  else
    output="${input}.sorted.desc.tsv.gz"
  fi
fi

tmpdir=${TMPDIR:-/tmp}
sort_mem=${SORT_MEM:-4G}

if [[ "$input" == *.gz ]]; then
  decompressor=(gzip -cd "$input")
else
  decompressor=(cat "$input")
fi

"${decompressor[@]}" \
  | {
      IFS= read -r header
      printf '%s\n' "$header"
      LC_ALL=C sort -T "$tmpdir" -S "$sort_mem" -t $'\t' -k4,4nr
    } \
  | gzip > "$output"

echo "$output"
