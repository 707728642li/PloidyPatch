#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from ploidypatch.wgdi_summary import parse_wgdi_collinearity


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collinearity", required=True, type=Path)
    parser.add_argument("--wgdi-gff", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shards", type=int, default=64)
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)
    (args.output_dir / "shards").mkdir()
    seqids: dict[str, str] = {}
    with args.wgdi_gff.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 6 or fields[1] in seqids:
                raise ValueError(f"Invalid target WGDI GFF line {line_number}")
            seqids[fields[1]] = fields[0]
    longest: dict[tuple[str, str], int] = defaultdict(int)
    for block in parse_wgdi_collinearity(args.collinearity, "target_self"):
        if len(block.pairs) < 20 or block.query_seqid == block.target_seqid:
            continue
        for left, right in block.pairs:
            if left == right or left not in seqids or right not in seqids:
                continue
            if seqids[left] == seqids[right]:
                continue
            pair = tuple(sorted((left, right)))
            longest[pair] = max(longest[pair], len(block.pairs))
    pairs = sorted(longest)
    if not pairs:
        raise ValueError("No target cross-chromosome structural pairs in blocks >=20")
    with (args.output_dir / "structural_pairs.tsv").open("x", encoding="utf-8") as handle:
        for left, right in pairs:
            handle.write(f"{left}\t{right}\n")
    shard_count = min(args.shards, len(pairs))
    for index in range(shard_count):
        with (args.output_dir / "shards" / f"pairs.{index:04d}.tsv").open(
            "x", encoding="utf-8"
        ) as handle:
            for left, right in pairs[index::shard_count]:
                handle.write(f"{left}\t{right}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
