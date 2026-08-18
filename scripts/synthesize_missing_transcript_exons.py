#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from ploidypatch.gff_compat import synthesize_missing_transcript_exons


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Append deterministic exon rows for transcripts lacking exons"
    )
    parser.add_argument("--input-gff", required=True)
    parser.add_argument("--output-gff", required=True)
    parser.add_argument(
        "--repair-parent-bounds",
        action="store_true",
        help=(
            "expand only gene/transcript bounds to same-seqid, same-strand "
            "child unions before synthesizing missing exons"
        ),
    )
    args = parser.parse_args()
    report = synthesize_missing_transcript_exons(
        args.input_gff,
        args.output_gff,
        repair_parent_bounds=args.repair_parent_bounds,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
