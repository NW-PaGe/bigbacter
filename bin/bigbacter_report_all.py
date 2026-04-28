#!/usr/bin/env python3
# bigbacter_report.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import glob
import os

from bigbacter_utils import logging_config

LOGGER = logging_config()


def combine_csvs(input_dir: str, timestamp: str):
    pattern = os.path.join(input_dir, "*.csv")
    files = glob.glob(pattern)

    if not files:
        raise ValueError(f"No CSV files found in directory: {input_dir}")

    all_columns = []
    all_rows = []

    for file in files:
        try:
            with open(file, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)

                # Track column order based on first appearance
                for col in (reader.fieldnames or []):
                    if col not in all_columns:
                        all_columns.append(col)

                all_rows.extend(rows)

        except Exception as e:
            print(f"Skipping file {file} due to error: {e}")

    if not all_rows and not all_columns:
        raise ValueError("No valid CSV files could be read.")

    output_file = f"{timestamp}-summary.csv"

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_columns, extrasaction="ignore")
        writer.writeheader()
        # Fill missing fields with empty string (mirrors reindex behavior)
        for row in all_rows:
            writer.writerow({col: row.get(col, "") for col in all_columns})

    print(f"Combined CSV written to: {output_file}")


def main():
    version = "1.0.0"
    parser = argparse.ArgumentParser(description="Combine CSV files into one.")
    parser.add_argument("input_dir", help="Directory containing CSV files")
    parser.add_argument("--timestamp", required=True, help="Timestamp for output filename")
    parser.add_argument("--version", action="version", version=version)

    args = parser.parse_args()

    LOGGER.info("%s v%s", os.path.basename(__file__).replace('.py', ''), version)
    LOGGER.info("Author: Jared Johnson")

    combine_csvs(args.input_dir, args.timestamp)


if __name__ == "__main__":
    main()