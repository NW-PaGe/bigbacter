#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Dict, Any


def open_maybe_gzip(path: str):
    """Open a text file that may be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def fasta_total_bases(ref_fasta: str) -> int:
    """Sum total bases across all FASTA records."""
    total = 0
    with open_maybe_gzip(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def extract_read_bases(fastp_json_path: str) -> int:
    """Extract summary.after_filtering.total_bases from fastp JSON."""
    with open_maybe_gzip(fastp_json_path) as f:
        data: Dict[str, Any] = json.load(f)

    try:
        return int(data["summary"]["after_filtering"]["total_bases"])
    except KeyError as e:
        raise ValueError(
            f"{fastp_json_path}: missing summary.after_filtering.total_bases"
        ) from e


def compute_rate(read_bases: int, ref_bases: int, max_depth: float) -> float:
    curr_depth = read_bases / ref_bases if ref_bases else 0
    if curr_depth <= 0:
        return 1.0
    return min(1.0, max_depth / curr_depth)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute downsampling rates from fastp JSON files."
    )
    p.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA (optionally gzipped).",
    )
    p.add_argument(
        "--max-depth",
        required=True,
        type=float,
        help="Target maximum depth.",
    )
    p.add_argument(
        "--out",
        default='ds_rate.csv',
        help="Output CSV file.",
    )
    p.add_argument(
        "fastp_json",
        nargs="+",
        help="One or more fastp JSON files (optionally gzipped).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    ref_bases = fasta_total_bases(args.ref)
    if ref_bases <= 0:
        raise SystemExit(f"{args.ref}: reference has 0 bases")

    rows = []
    for jp in args.fastp_json:
        sid = Path(jp).name   # 👈 simple, explicit, predictable
        read_bases = extract_read_bases(jp)
        rate = compute_rate(read_bases, ref_bases, args.max_depth)

        rows.append({
            "filename": sid,
            "ref_bases": ref_bases,
            "read_bases": read_bases,
            "rate": f"{rate:.6f}",
        })

    rows.sort(key=lambda r: r["filename"])

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["filename", "ref_bases", "read_bases", "rate"]
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
