from __future__ import annotations

import argparse
import json

from ploidypatch.gff_compat import synthesize_root_transcript_genes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gff", required=True)
    parser.add_argument("--output-gff", required=True)
    parser.add_argument("--mapping-tsv", required=True)
    args = parser.parse_args()
    report = synthesize_root_transcript_genes(
        args.input_gff, args.output_gff, args.mapping_tsv
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
