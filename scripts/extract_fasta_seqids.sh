#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 INPUT_FASTA[.gz] OUTPUT_FASTA SEQID_CSV" >&2
  exit 2
fi

input_fasta=$1
output_fasta=$2
seqid_csv=$3
temporary_output="${output_fasta}.tmp.$$"

cleanup() {
  rm -f -- "$temporary_output"
}
trap cleanup EXIT

case "$input_fasta" in
  *.gz) reader=(gzip -dc -- "$input_fasta") ;;
  *) reader=(cat -- "$input_fasta") ;;
esac

"${reader[@]}" | awk -v wanted_csv="$seqid_csv" '
BEGIN {
  count = split(wanted_csv, seqids, ",")
  for (i = 1; i <= count; i++) {
    wanted[seqids[i]] = 1
  }
}
/^>/ {
  seqid = substr($1, 2)
  keep = (seqid in wanted)
}
keep { print }
' > "$temporary_output"

observed=$(grep -c '^>' "$temporary_output")
expected=$(awk -F, '{ print NF }' <<< "$seqid_csv")
if [[ "$observed" -ne "$expected" ]]; then
  echo "expected $expected FASTA records but extracted $observed" >&2
  exit 1
fi

mv -- "$temporary_output" "$output_fasta"
trap - EXIT
