#!/usr/bin/env python3
# bigbacter_report_cluster.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov


from typing import List, Dict, Any, Tuple, Optional

import json
import csv
import argparse
import logging
from typing import List, Dict, Any
from Bio import Phylo
import numpy as np
from sklearn.cluster import DBSCAN
from pathlib import Path
from datetime import datetime, timezone
import gzip
import os

from bigbacter_utils import logging_config, get_assembly_stem

LOGGER = logging_config()

COL_ORDER = [
    'sample', 'run', 'status', 'included', 
    'taxa', 'cluster', 'partition', 'strong_links', 
    'inter_links', 'genome_fraction', 'core_fraction', 
    'length', 'masked', 'missing', 'mixed', 'variants', 
    'recomb_masked'
]

# -------------------------------
#  HELPERS
# -------------------------------
def _rename_reference(data: Dict[str, Dict[str, Any]], ref_name: str, depth: int = 1) -> Dict[str, Dict[str, Any]]:
    """Rename 'Reference' row and column keys to ref_name.
    
    depth=1: rename outer row keys only
    depth=2: rename outer row keys and inner column keys (default)
    """
    if ref_name == "Reference":
        return data
    if depth == 2:
        for v in data.values():
            if isinstance(v, dict) and "Reference" in v:
                v[ref_name] = v.pop("Reference")
    if "Reference" in data:
        data[ref_name] = data.pop("Reference")

def _ensure_path(p: Path, label: str) -> bool:
    if not p.exists():
        logging.warning(f"{label} missing: {p}")
        return False
    return True

def _attach_text_file(mr: Dict[str, Any], id: str, path: Path) -> None:
    mr.setdefault('files', {}).setdefault(id, {})
    mr['files'][id]['blob'] = path.read_text(encoding='utf-8')
    mr['files'][id]['name'] = path.name

def _delete_pane_by_id(node, pane_id):
    """
    Recursively search through 'children' lists and delete
    the dict whose 'id' matches pane_id.

    Returns True if deleted, False otherwise.
    """
    if not isinstance(node, dict):
        return False

    children = node.get("children", [])
    for i, child in enumerate(children):
        if child.get("id") == pane_id:
            del children[i]
            return True

        if _delete_pane_by_id(child, pane_id):
            return True

    return False

def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def _extend_dict(main: Dict[str, Dict[str, Any]], new: Dict[str, Dict[str, Any]]) -> None:
    """Merge nested dict-of-dicts by key."""
    for k, v in new.items():
        main[k] = (main.get(k, {}) | v)

def _list2map(data: List[Dict[str, Any]], key: str) -> Dict[str, Dict[str, Any]]:
    return {
        row[key]: {k: v for k, v in row.items() if k != key}
        for row in data
        if key in row and row[key] not in (None, "")
    }

def _sniff_open(path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', encoding='utf-8-sig')
    else:
        return open(path, 'r', encoding='utf-8-sig')
    
def _sniff_delimiter(path):
    # Robust delimiter sniff with fallback
    sample = path.read(2048)
    path.seek(0)
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(sample)
        return dialect.delimiter
    except Exception:
        return ','
    
def _load_delim(path: Path) -> List[Dict[str, str]]:
    with _sniff_open(path) as f:
        return list(csv.DictReader(f, delimiter=_sniff_delimiter(f)))    

def _find_links(data, strong_thresh: float, inter_thresh: float) -> Dict[str, Dict[str, Any]]:
    res = {}
    for k1, v1 in data.items():
        strong_links, inter_links = [], []
        for k2, v2 in v1.items():
            if float(v2) <= float(strong_thresh):
                strong_links.append(k2)
            elif float(v2) <= float(inter_thresh):
                inter_links.append(k2)
        res[k1] = {'inter_links': inter_links, 'strong_links': strong_links}
    return res

def _patristic_distance_matrix(tree) -> Tuple[List[str], np.ndarray]:
    tips = tree.get_terminals()
    labels = [t.name for t in tips]
    n = len(tips)
    M = np.zeros((n, n), dtype=float)
    for i, a in enumerate(tips):
        for j in range(i, n):
            d = tree.distance(a, tips[j])
            M[i, j] = M[j, i] = d
    return labels, M

def _partition_tree(D: np.ndarray, eps: float, min_samples: int = 1) -> List[int]:
    db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    y = db.fit_predict(D)

    # Keep your "clusters are 1-based, and noise becomes unique clusters" behavior
    max_cluster = max(y) if len(y) else -1
    clusters = [c + 1 for c in y]
    for i, c in enumerate(clusters):
        if c == 0:
            c = max_cluster + 2
            max_cluster = c - 1
        clusters[i] = c
    return clusters

def _dist2csv(data, outfile):
    # Collect all second-level keys
    samples = sorted(data.keys())

    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[''] + samples)
        writer.writeheader()

        for s in samples:

            if s not in data:
                continue

            row = data[s].copy()
            if s not in row:
                row[s] = 0
            row[''] = s
        
            writer.writerow(row)

def _find_file(indir: Path, pattern: str) -> Optional[Path]:
    """Find first file matching glob pattern in directory."""
    matches = list(indir.glob(pattern))
    if matches:
        LOGGER.info(f"Found {pattern}: {matches[0]}")
        return matches[0]
    LOGGER.warning(f"No file matching pattern '{pattern}' in {indir}")
    return None

# -------------------------------
#  CORE LOGIC
# -------------------------------

def process_run(indir: Path, args: argparse.Namespace) -> None:
    out: Dict[str, Dict[str, Any]] = {}

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    try:
        timestamp_epoch = int(args.timestamp)
        timestamp_readable = datetime.fromtimestamp(timestamp_epoch, tz=timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    except (ValueError, TypeError):
        # If timestamp is not a valid epoch, use as-is
        timestamp_epoch = args.timestamp
        timestamp_readable = None

    # ----------------------------
    # Find input files by pattern
    # ----------------------------
    samplesheet_file    = indir / "samplesheet.csv"
    snp_summary_file    = _find_file(indir, "summary.csv")
    snp_dist_file       = _find_file(indir, "dist_wide.csv")
    snp_tree_file       = _find_file(indir, "core.aln.*")
    mh_dist_file        = _find_file(indir, "*.tsv.gz")
    ref_meta_file       = indir / "ref.json"
    microreact_template = indir / "template.microreact"

    # ----------------------------
    # Load samplesheet
    # ----------------------------
    samples = []
    if samplesheet_file.exists():
        samplesheet = _load_delim(samplesheet_file)
        samples = [rec['sample'] for rec in samplesheet if rec.get('sample')]
        LOGGER.info(f"Loaded samplesheet with {len(samples)} samples")
    else:
        LOGGER.warning(f"Samplesheet not found: {samplesheet_file}")

    # ----------------------------
    # Reference metadata
    # ----------------------------
    ref_name = "Reference"
    if ref_meta_file:
        try:
            ref_meta = load_json(ref_meta_file)
            ref_origin = get_assembly_stem(ref_meta.get('name', '')).strip()
            if ref_origin:
                ref_name = f"{ref_name}_[{ref_origin}]"
        except:
            LOGGER.warning("Failed to load reference metadata")

    # ----------------------------
    # Alignment stats
    # ----------------------------
    stats_dict: Dict[str, Dict[str, Any]] = {}
    if snp_summary_file:
        stats_dict = _list2map(_load_delim(snp_summary_file), key='name')
        if stats_dict:
            _rename_reference(stats_dict, ref_name)
            _extend_dict(out, stats_dict)
        LOGGER.info(f"Loaded alignment stats ({len(stats_dict)} rows) from {snp_summary_file}")
    else:
        LOGGER.warning("No summary.csv found")

    # ----------------------------
    # SNP distance matrix
    # ----------------------------
    snp_dist_out = Path(outdir) / f"{args.prefix}_snp-dist{args.suffix}csv"
    if snp_dist_file:
        snp_dist_data = _list2map(_load_delim(snp_dist_file), key='name')
        _rename_reference(snp_dist_data, ref_name, depth=2)
        _dist2csv(snp_dist_data, outfile=snp_dist_out)
        _extend_dict(out, _find_links(snp_dist_data, args.strong_link, args.inter_link))
        LOGGER.info(f"Processed SNP distance matrix from {snp_dist_file}")

    # ----------------------------
    # MinHash containment distance
    # ----------------------------
    mh_dist_out = Path(outdir) / f"{args.prefix}_db-dist{args.suffix}csv"
    if mh_dist_file:
        mh_dist_data = _list2map(_load_delim(mh_dist_file), key='')
        for k1, v1 in mh_dist_data.items():
            for k2 in v1:
                mh_dist_data[k1][k2] = round(100 * float(v1[k2]), 2)

        _dist2csv(mh_dist_data, outfile=mh_dist_out)
        LOGGER.info(f"Processed MinHash containment distance matrix from {mh_dist_file}")


    # ----------------------------
    # Tree + partitions
    # ----------------------------
    snp_tree_out = Path(outdir) / f"{args.prefix}_snp-tree{args.suffix}nwk"
    if snp_tree_file:
        tree = Phylo.read(snp_tree_file, 'newick')
        tree.root_at_midpoint()

        # Attempt to scale by reference length if present
        scale = 1.0
        try:
            vb = stats_dict.get(ref_name, {}).get('length')
            if vb not in (None, ""):
                scale = float(vb)
            else:
                LOGGER.warning("Branch lengths not scaled (Reference length missing)")
        except Exception:
            LOGGER.warning("Branch lengths not scaled (failed to parse Reference length)")

        for clade in tree.find_clades():
            # Update reference name
            if clade.is_terminal() and clade.name == 'Reference':
                clade.name = ref_name
            # Rescale branch
            bl = getattr(clade, "branch_length", None)
            if bl is not None:
                try:
                    clade.branch_length = float(bl) * scale
                except Exception:
                    pass

        Phylo.write(tree, snp_tree_out, "newick")

        labels, dist = _patristic_distance_matrix(tree)

        parts = _partition_tree(dist, eps=args.partition_distance)
        _extend_dict(out, {l: {'partition': p} for l, p in zip(labels, parts)})

        LOGGER.info("Processed tree with %d tips and partitions; wrote %s", len(labels), snp_tree_out)
            
    # ----------------------------
    # Build final summary and write CSV
    # ----------------------------
    summary_out = Path(outdir) / f"{args.prefix}_summary{args.suffix}csv"
    summary_records: List[Dict[str, Any]] = []

    for sample_id, sample_data in out.items():
        rec: Dict[str, Any] = {
            "sample": ref_name if sample_id == "Reference" else sample_id,
            "status": "new" if sample_id in samples else "old",
            "run": timestamp_epoch,
            "taxa": args.taxa,
            "cluster": args.cluster,
            "recomb_masked": args.recomb_masked
        }
        for k, v in sample_data.items():
            if k in {"name"}:
                continue
            rec[k] = ":".join(map(str, v)) if isinstance(v, list) else v
        summary_records.append(rec)

    if not summary_records:
        LOGGER.warning("No records to write")
        return

    # ---- column ordering ----
    all_keys: set[str] = set()
    for rec in summary_records:
        all_keys.update(rec.keys())

    # 1) take COL_ORDER in-order (only those present)
    ordered = [c for c in COL_ORDER if c in all_keys]

    # 2) add any remaining columns not specified (sorted for stability)
    remaining = sorted(all_keys - set(ordered))

    summary_fieldnames = ordered + remaining

    with open(summary_out, "w", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=summary_fieldnames,
            delimiter=","
        )
        writer.writeheader()
        writer.writerows(summary_records)

    LOGGER.info(
        f"Wrote {summary_out} with {len(summary_records)} records and {len(summary_fieldnames)} columns"
    )

    # ----------------------------
    # Microreact packaging (optional)
    # ----------------------------
    if microreact_template.exists():

        mr_map = {
            "summary_file": {"pane_id": "table-1", "file": summary_out},
            "snp_tree_file": {"pane_id": "tree-1", "file": snp_tree_out},
            "snp_dist_file": {"pane_id": "matrix-1", "file": snp_dist_out},
            "mh_dist_file": {"pane_id": "matrix-2", "file": mh_dist_out},
        }

        mr_json = load_json(microreact_template)

        layout = mr_json['panes']['model']['layout']

        for k, v in mr_map.items():
            pane_id, filepath = v.values()
            if _ensure_path(filepath, k):
                _attach_text_file(mr_json, k, filepath)
            else:
                _delete_pane_by_id(layout, pane_id)

        # Meta
        # Build description
        desc_parts = [
            "BigBacter Microreact Report",
            f"Run time: {timestamp_readable} ({timestamp_epoch})",
            f"Species: {args.taxa}",
            f"Cluster: {args.cluster}",
            f"Recombination masked: {args.recomb_masked}"
        ]
        
        mr_json['meta']['description'] = " | ".join(desc_parts)

        name_parts = [args.taxa, args.cluster, "[masked]" if args.recomb_masked else "", f"({timestamp_readable})"]
        mr_json['meta']['name'] = ' '.join(name_parts)
        mr_json['meta']['timestamp'] = (
            datetime.now(timezone.utc)
            .isoformat(timespec="milliseconds")
            .replace("+00:00", "Z")
        )

        # Columns from written summary (if exists)
        if summary_fieldnames:
            mr_json['tables']['table-1']['columns'] = [{"field": f, "fixed": False} for f in summary_fieldnames]

        microreact_out = Path(outdir) / f"{args.prefix}{args.suffix}microreact"
        with open(microreact_out, "w", encoding="utf-8") as f:
            json.dump(mr_json, f, indent=2, ensure_ascii=False)
        LOGGER.info("Wrote %s", microreact_out)
    else:
        LOGGER.warning(f"Microreact template not found: {microreact_template}")

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    version = "1.0.0"
    parser = argparse.ArgumentParser(
        description="Summarize BigBacter run outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected input directory structure:
  - samplesheet.csv (required)
  - summary.csv (SNP summary stats)
  - dist_wide.csv (SNP distance matrix)
  - core.aln.* (Newick tree file)
  - *.tsv.gz (K-mer distance matrix)
  - template.microreact (optional, for Microreact output)
        """
    )
    parser.add_argument("indir", type=str, help="Input directory containing all required files")
    parser.add_argument("--timestamp", type=str, required=True, help="Run timestamp (epoch)")
    parser.add_argument("--taxa", type=str, required=True, help="Taxonomy/species")
    parser.add_argument("--cluster", type=str, required=True, help="Cluster identifier")
    parser.add_argument("--recomb-masked", action="store_true", help="Recombination masking was performed")
    parser.add_argument("--partition-distance", default=100.0, type=float, help="Distance threshold for tree partitioning")
    parser.add_argument("--strong-link", default=10.0, type=float, help="SNP distance threshold for strong links")
    parser.add_argument("--inter-link", default=50.0, type=float, help="SNP distance threshold for intermediate links")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory")
    parser.add_argument("--version", action="version", version=version)
    args = parser.parse_args()

    LOGGER.info("%s v%s", os.path.basename(__file__).replace('.py', ''), version)
    LOGGER.info("Author: Jared Johnson")

    # Build prefix from components
    prefix_parts = [args.timestamp, args.taxa, args.cluster]
    suffix = ".masked." if args.recomb_masked else "."

    args.prefix = "-".join(prefix_parts)
    args.suffix = suffix
    
    LOGGER.info(f"Output prefix + suffix: {args.prefix + args.suffix}")

    indir = Path(args.indir)
    if not indir.exists():
        LOGGER.error(f"Input directory does not exist: {indir}")
        raise FileNotFoundError(f"Input directory not found: {indir}")

    if not indir.is_dir():
        LOGGER.error(f"Input path is not a directory: {indir}")
        raise NotADirectoryError(f"Not a directory: {indir}")

    process_run(indir, args)

if __name__ == "__main__":
    main()