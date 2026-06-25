#!/usr/bin/env python3

import csv

MATRIX = "dist_wide.csv"
SUMMARY = "summary.csv"

stats = {}
with open(MATRIX) as f:
    reader = csv.reader(f)
    header = next(reader)
    col_names = header[1:]  # sample names from header (excluding first cell)

    for row_idx, row in enumerate(reader):
        name = row[0]
        if name == "Reference":   # don't produce stats for the reference itself
            continue
        values = []
        for i in range(1, len(row)):
            col_name = col_names[i - 1]
            if col_name == name:          # skip diagonal
                continue
            if col_name == "Reference":   # skip reference column
                continue
            try:
                values.append(int(row[i]))
            except ValueError:
                continue
        if values:
            stats[name] = {
                "MEAN_SNP_DIST": round(sum(values) / len(values)),
                "MIN_SNP_DIST": min(values),
                "MAX_SNP_DIST": max(values),
            }

# Merge into summary.csv
stat_cols = ["MEAN_SNP_DIST", "MIN_SNP_DIST", "MAX_SNP_DIST"]

try:
    with open(SUMMARY) as f:
        reader = csv.DictReader(f)
        existing_fields = reader.fieldnames or []
        rows = list(reader)

    new_fields = existing_fields + [c for c in stat_cols if c not in existing_fields]
    for row in rows:
        for col in stat_cols:
            row[col] = ""  # clear old values before updating
        if row["name"] in stats:
            row.update(stats[row["name"]])

    with open(SUMMARY, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=new_fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done — updated {SUMMARY}")

except FileNotFoundError:
    with open(SUMMARY, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["name"] + stat_cols)
        writer.writeheader()
        for name, s in stats.items():
            writer.writerow({"name": name, **s})

    print(f"No existing {SUMMARY} found — created it with SNP stats only.")