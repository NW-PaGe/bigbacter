#!/usr/bin/env python3
# bigbacter_report.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import glob
import os

from bigbacter_utils import logging_config

LOGGER = logging_config()

MASK_COLUMN = "recomb_masked"
MASK_VALUES = {"True": True, "False": False}


def load_csvs(input_dir: str):
    pattern = os.path.join(input_dir, "*.csv")
    files = sorted(glob.glob(pattern))

    if not files:
        raise ValueError(f"No CSV files found in directory: {input_dir}")

    all_columns = []
    all_rows = []

    for file in files:
        try:
            with open(file, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames or []
                rows = list(reader)
        except Exception as e:
            LOGGER.warning("Skipping file %s due to error: %s", file, e)
            continue

        if MASK_COLUMN not in fieldnames:
            raise ValueError(f"Required column '{MASK_COLUMN}' not found in file: {file}")

        # Track column order based on first appearance
        for col in fieldnames:
            if col not in all_columns:
                all_columns.append(col)

        # Tag each row with its source file for error reporting
        for line_num, row in enumerate(rows, start=2):
            all_rows.append((file, line_num, row))

    if not all_rows and not all_columns:
        raise ValueError("No valid CSV files could be read.")

    return all_columns, all_rows


def split_rows(rows):
    masked, unmasked = [], []

    for file, line_num, row in rows:
        value = row.get(MASK_COLUMN)
        if value not in MASK_VALUES:
            raise ValueError(
                f"Invalid '{MASK_COLUMN}' value {value!r} in {file} (line {line_num}); "
                f"expected exactly 'True' or 'False'."
            )
        (masked if MASK_VALUES[value] else unmasked).append(row)

    return masked, unmasked


def write_csv(output_file: str, columns, rows):
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        # Fill missing fields with empty string (mirrors reindex behavior)
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})
    LOGGER.info("Wrote %d row(s) to: %s", len(rows), output_file)


def combine_csvs(input_dir: str, timestamp: str):
    columns, rows = load_csvs(input_dir)
    masked, unmasked = split_rows(rows)

    # Header-only files are still written so downstream steps always find both outputs
    write_csv(f"{timestamp}-summary.csv", columns, unmasked)
    write_csv(f"{timestamp}-summary.masked.csv", columns, masked)


def main():
    version = "1.1.0"
    parser = argparse.ArgumentParser(
        description="Combine CSV files into masked and unmasked summaries based on the 'recomb_masked' column."
    )
    parser.add_argument("input_dir", help="Directory containing CSV files")
    parser.add_argument("--timestamp", required=True, help="Timestamp for output filenames")
    parser.add_argument("--version", action="version", version=version)

    args = parser.parse_args()

    LOGGER.info("%s v%s", os.path.basename(__file__).replace('.py', ''), version)
    LOGGER.info("Author: Jared Johnson")

    combine_csvs(args.input_dir, args.timestamp)


if __name__ == "__main__":
    main()