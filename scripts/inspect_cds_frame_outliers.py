#!/usr/bin/env python3
"""Print representative CDS records whose sequence length is not divisible by 3."""

from __future__ import annotations

import argparse

from ploidypatch.io import fasta_relation_id, iter_fasta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cds_fasta")
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()

    emitted = 0
    for primary_id, header, sequence in iter_fasta(args.cds_fasta):
        if len(sequence) % 3 == 0:
            continue
        transcript_id, source = fasta_relation_id(primary_id, header)
        print(
            transcript_id,
            primary_id,
            len(sequence),
            len(sequence) % 3,
            source,
            header,
            sep="\t",
        )
        emitted += 1
        if emitted >= args.limit:
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
