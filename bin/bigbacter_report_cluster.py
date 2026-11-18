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
    'id', 'run', 'status', 'included', 
    'taxa', 'cluster', 'partition',
    'accessory_partition', 'combined_partition',
    'strong_links',
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
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
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

def _get_delimiter(path: Path) -> str:
    """Determine delimiter from file extension."""
    suffixes = path.suffixes  # e.g. ['.tsv', '.gz'] or ['.csv']
    for s in suffixes:
        if s in ('.tsv', '.txt'):
            return '\t'
        if s == '.csv':
            return ','
    return ','  # default fallback

def _sniff_open(path: Path):
    if '.gz' in path.suffixes:
        return gzip.open(path, 'rt', encoding='utf-8-sig')
    else:
        return open(path, 'r', encoding='utf-8-sig')

def _load_delim(path: Path) -> List[Dict[str, str]]:
    delimiter = _get_delimiter(path)
    with _sniff_open(path) as f:
        return list(csv.DictReader(f, delimiter=delimiter))

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

def _partition_distances(D: np.ndarray, eps: float, min_samples: int = 1) -> List[int]:
    """Single-linkage-style grouping of a precomputed distance matrix via DBSCAN.

    Used for both the core (SNP patristic) and accessory (event count) matrices.
    eps is in the units of D, so the two callers pass different thresholds.
    """
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
    samples = sorted(data.keys())

    # Collect ALL column keys across all rows (inner keys may differ from outer)
    all_cols = sorted({col for row in data.values() for col in row.keys()})

    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[''] + all_cols, extrasaction='ignore')
        writer.writeheader()

        for s in samples:
            row = data[s].copy()
            row[s] = row.get(s, 0)   # self-distance = 0
            row[''] = s
            writer.writerow(row)

def _dist2csv_long(data, outfile, value_col: str = 'dist') -> None:
    """Write a dict-of-dicts distance matrix in long (tidy) format.

    Output columns: id1, id2, <value_col>
    Self-distances (diagonal) are emitted as 0 when not present in the data.
    """
    samples = sorted(data.keys())
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(['id1', 'id2', value_col])
        for s1 in samples:
            row = data[s1]
            for s2 in sorted(set(row.keys()) | {s1}):
                val = row.get(s2, 0) if s2 == s1 else row.get(s2)
                writer.writerow([s1, s2, val])

def _load_rtab(path: Path) -> Tuple[List[str], np.ndarray]:
    """Load a Panaroo/Roary gene presence-absence Rtab file.

    Returns the sample names and a samples x genes 0/1 array.
    """
    with _sniff_open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        samples = next(reader)[1:]
        rows = [[1 if v.strip() not in ('', '0') else 0 for v in row[1:]] for row in reader if row]
    if any(len(r) != len(samples) for r in rows):
        raise ValueError(f"{path}: rows do not all have {len(samples)} sample columns")
    pa = np.array(rows, dtype=np.uint8).reshape(len(rows), len(samples))
    return samples, pa.T

def _gene_distance_matrix(pa: np.ndarray) -> np.ndarray:
    """Pairwise number of genes present in one sample but not the other (samples x genes 0/1 input)."""
    x = pa.astype(np.float64)
    shared = x @ x.T
    counts = np.diag(shared)
    return np.rint(counts[:, None] + counts[None, :] - 2 * shared).astype(int)

def _load_square_matrix(path: Path) -> Tuple[List[str], np.ndarray]:
    """Load a square wide-format distance matrix (first column is the sample labels).

    Used for the pan-dists accessory event distances. Returns labels and an
    integer matrix in the order the labels appear.
    """
    delimiter = _get_delimiter(path)
    with _sniff_open(path) as f:
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)[1:]
        labels, rows = [], []
        for row in reader:
            if not row:
                continue
            labels.append(row[0].strip())
            rows.append([float(v) if v.strip() not in ('', 'NA') else 0.0 for v in row[1:]])
    header = [h.strip() for h in header]
    if labels != header:
        raise ValueError(f"{path}: row labels and column headers differ")
    if any(len(r) != len(labels) for r in rows):
        raise ValueError(f"{path}: rows do not all have {len(labels)} columns")
    return labels, np.rint(np.array(rows, dtype=float)).astype(int)

def _neighbor_joining(labels: List[str], dist: np.ndarray):
    """Build a neighbor-joining tree (Saitou & Nei 1987) from a distance matrix.

    Negative branch lengths are set to zero and the difference is moved to the sister branch.
    Returns an unrooted Bio.Phylo tree. Needs at least two samples.
    """
    from Bio.Phylo.BaseTree import Clade, Tree

    if len(labels) < 2:
        raise ValueError("Neighbor joining needs at least two samples")

    D = np.array(dist, dtype=float)
    nodes = [Clade(name=str(l)) for l in labels]

    if len(nodes) == 2:
        nodes[0].branch_length = nodes[1].branch_length = max(D[0, 1], 0.0) / 2
        return Tree(root=Clade(clades=nodes), rooted=False)

    while len(nodes) > 3:
        m = len(nodes)
        r = D.sum(axis=1)
        Q = (m - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(Q.argmin()), m)
        if i > j:
            i, j = j, i

        dij = D[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2))
        lj = dij - li
        if li < 0:
            li, lj = 0.0, max(lj - li, 0.0)
        elif lj < 0:
            li, lj = max(li - lj, 0.0), 0.0

        nodes[i].branch_length = li
        nodes[j].branch_length = lj
        parent = Clade(clades=[nodes[i], nodes[j]])

        # distances from the new node replace row/column i; row/column j is dropped
        dk = 0.5 * (D[i] + D[j] - dij)
        D[i, :] = dk
        D[:, i] = dk
        D[i, i] = 0.0
        D = np.delete(np.delete(D, j, axis=0), j, axis=1)
        nodes[i] = parent
        del nodes[j]

    # the last three nodes meet at the root
    nodes[0].branch_length = max((D[0, 1] + D[0, 2] - D[1, 2]) / 2, 0.0)
    nodes[1].branch_length = max((D[0, 1] + D[1, 2] - D[0, 2]) / 2, 0.0)
    nodes[2].branch_length = max((D[0, 2] + D[1, 2] - D[0, 1]) / 2, 0.0)
    return Tree(root=Clade(clades=nodes), rooted=False)

def _root_at_midpoint(tree):
    """Root a tree halfway along the path between its two most distant tips.

    Returns a new Bio.Phylo tree with the same tip names and tip-to-tip distances. Used instead of
    Bio.Phylo's root_at_midpoint, which lengthens the tree when the midpoint is on a tip's branch.
    """
    from Bio.Phylo.BaseTree import Clade, Tree

    # the tree as an undirected graph
    adj: Dict[Any, List[Tuple[Any, float]]] = {}
    for parent in tree.find_clades(order='level'):
        adj.setdefault(parent, [])
        for child in parent.clades:
            length = float(child.branch_length or 0.0)
            adj[parent].append((child, length))
            adj.setdefault(child, []).append((parent, length))
    tips = [c for c in adj if c.is_terminal()]

    def _sweep(start):
        dist, prev, stack = {start: 0.0}, {start: None}, [start]
        while stack:
            u = stack.pop()
            for v, length in adj[u]:
                if v not in dist:
                    dist[v], prev[v] = dist[u] + length, u
                    stack.append(v)
        return dist, prev

    # the two most distant tips (a, b), then the branch (u, v) on their path that holds the midpoint
    dist, _ = _sweep(tips[0])
    a = max(tips, key=lambda t: dist[t])
    dist, prev = _sweep(a)
    b = max((t for t in tips if t is not a), key=lambda t: dist[t])
    half = dist[b] / 2
    v = b
    while dist[prev[v]] > half:
        v = prev[v]
    u = prev[v]

    def _build(node, came_from, length):
        top = Clade(name=node.name if node.is_terminal() else None, branch_length=length)
        stack = [(node, came_from, top)]
        while stack:
            n, source, new = stack.pop()
            for nb, nb_length in adj[n]:
                if nb is source:
                    continue
                child = Clade(name=nb.name if nb.is_terminal() else None, branch_length=nb_length)
                new.clades.append(child)
                stack.append((nb, n, child))
        return top

    root = Clade(clades=[_build(u, v, half - dist[u]), _build(v, u, dist[v] - half)])

    # merge nodes left with a single child (the old root, when it had two children)
    stack = [root]
    while stack:
        clade = stack.pop()
        for idx, child in enumerate(clade.clades):
            while len(child.clades) == 1:
                only = child.clades[0]
                only.branch_length = (only.branch_length or 0.0) + (child.branch_length or 0.0)
                clade.clades[idx] = child = only
            stack.append(child)

    return Tree(root=root, rooted=True)

def _write_newick(tree, outfile: Path) -> None:
    """Write Newick without Biopython's automatic single-quote wrapping."""
    from io import StringIO
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    with open(outfile, "w") as fh:
        fh.write(buf.getvalue().replace("'", ""))

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
    gene_pa_file        = _find_file(indir, "*.Rtab")
    acc_dist_file       = _find_file(indir, "dist_events_all.tsv")
    ref_meta_file       = next((f for f in (indir / "ref.json", indir / "ref.json.gz") if f.exists()), None)
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
                ref_name = f"{ref_name}_{ref_origin}"
        except Exception as e:
            LOGGER.warning(f"Failed to load reference metadata from {ref_meta_file}: {e}")
    else:
        LOGGER.warning("No reference metadata (ref.json) found")

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
    snp_dist_out = Path(outdir) / f"{args.prefix}_snp-dist-wide{args.suffix}csv"
    snp_dist_long_out = Path(outdir) / f"{args.prefix}_snp-dist-long{args.suffix}csv"
    if snp_dist_file:
        snp_dist_data = _list2map(_load_delim(snp_dist_file), key='name')
        _rename_reference(snp_dist_data, ref_name, depth=2)
        _dist2csv(snp_dist_data, outfile=snp_dist_out)
        _dist2csv_long(snp_dist_data, outfile=snp_dist_long_out)
        _extend_dict(out, _find_links(snp_dist_data, args.strong_link, args.inter_link))
        LOGGER.info(f"Processed SNP distance matrix from {snp_dist_file}")

    # ----------------------------
    # MinHash containment distance
    # ----------------------------
    mh_dist_out = Path(outdir) / f"{args.prefix}_mh-dist-wide{args.suffix}csv"
    mh_dist_long_out = Path(outdir) / f"{args.prefix}_mh-dist-long{args.suffix}csv"
    if mh_dist_file:
        mh_dist_data = _list2map(_load_delim(mh_dist_file), key='')
        LOGGER.info(f"mh_dist_data keys (first 5): {list(mh_dist_data.keys())[:5]}")
        LOGGER.info(f"out keys (first 5): {list(out.keys())[:5]}")
        LOGGER.info(f"mh_dist_data sample inner keys (first row): {list(list(mh_dist_data.values())[0].keys())[:5] if mh_dist_data else 'EMPTY'}")
        for k1, v1 in mh_dist_data.items():
            for k2 in v1:
                mh_dist_data[k1][k2] = round(100 * float(v1[k2]), 2)

        # Filter to only samples that will appear in the summary
        mh_dist_data = {
            k1: {k2: v for k2, v in v1.items() if k2 in out}
            for k1, v1 in mh_dist_data.items()
            if k1 in out
        }

        _dist2csv(mh_dist_data, outfile=mh_dist_out)
        _dist2csv_long(mh_dist_data, outfile=mh_dist_long_out, value_col='distance')
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

        _write_newick(tree, snp_tree_out)

        labels, dist = _patristic_distance_matrix(tree)

        parts = _partition_distances(dist, eps=args.partition_distance)
        _extend_dict(out, {l: {'partition': p} for l, p in zip(labels, parts)})

        LOGGER.info("Processed tree with %d tips and partitions; wrote %s", len(labels), snp_tree_out)
            
    # ----------------------------
    # Pairwise accessory differences + neighbor-joining tree (optional)
    #
    # Preferred source is pan-dists dist_events_all.tsv: accessory genes lying within
    # N genes of one another are grouped into gain/loss events, and each pair is scored
    # by how many events differ (multi-cluster blocks + isolated single genes). A single
    # plasmid acquisition therefore counts once rather than once per gene it carries.
    #
    # Falls back to counting differing genes straight from the Panaroo Rtab when
    # pan-dists output is not present, so runs without the PAN_DISTS process still work.
    # ----------------------------
    acc_dist_out = Path(outdir) / f"{args.prefix}_accessory-dist-wide{args.suffix}csv"
    acc_dist_long_out = Path(outdir) / f"{args.prefix}_accessory-dist-long{args.suffix}csv"
    acc_tree_out = Path(outdir) / f"{args.prefix}_accessory-tree{args.suffix}nwk"

    acc_samples: List[str] = []
    acc_dist: Optional[np.ndarray] = None
    acc_source = None

    if acc_dist_file:
        try:
            acc_samples, acc_dist = _load_square_matrix(acc_dist_file)
            acc_source = "accessory events"
            LOGGER.info(
                "Processed accessory event distances (%d samples) from %s",
                len(acc_samples), acc_dist_file
            )
        except Exception as e:
            LOGGER.warning("Could not read %s (%s); falling back to gene differences", acc_dist_file, e)
            acc_dist = None

    if acc_dist is None and gene_pa_file:
        acc_samples, gene_pa = _load_rtab(gene_pa_file)
        acc_dist = _gene_distance_matrix(gene_pa)
        acc_source = "gene differences"
        LOGGER.info(
            "Processed gene presence-absence matrix (%d samples, %d genes) from %s",
            len(acc_samples), gene_pa.shape[1], gene_pa_file
        )

    if acc_dist is not None:
        acc_dist_data = {
            s1: {s2: int(acc_dist[i, j]) for j, s2 in enumerate(acc_samples)}
            for i, s1 in enumerate(acc_samples)
        }
        _dist2csv(acc_dist_data, outfile=acc_dist_out)
        _dist2csv_long(acc_dist_data, outfile=acc_dist_long_out)
        LOGGER.info("Pairwise accessory differences measured as %s; wrote %s", acc_source, acc_dist_out)

        # Partition on the event distances themselves rather than on patristic
        # distances from the NJ tree: the matrix is the primary object here, so
        # going through a tree reconstruction would only add distortion.
        if acc_source == "gene differences":
            LOGGER.warning(
                "Accessory partitions are being computed from gene differences, not "
                "events; --acc-partition-distance %s is in event units and is very "
                "likely wrong for this input",
                args.acc_partition_distance
            )
        acc_parts = _partition_distances(
            acc_dist.astype(float), eps=args.acc_partition_distance
        )
        _extend_dict(out, {l: {'accessory_partition': p} for l, p in zip(acc_samples, acc_parts)})
        LOGGER.info(
            "Assigned %d accessory partitions across %d samples (threshold %s %s)",
            len(set(acc_parts)), len(acc_samples), args.acc_partition_distance, acc_source
        )

        if len(acc_samples) > 1:
            # sorted so the tree does not depend on the column order of the input
            order = sorted(range(len(acc_samples)), key=lambda i: acc_samples[i])
            acc_tree = _neighbor_joining(
                [acc_samples[i] for i in order],
                acc_dist[np.ix_(order, order)]
            )
            _write_newick(_root_at_midpoint(acc_tree), acc_tree_out)
            LOGGER.info("Built accessory neighbor-joining tree; wrote %s", acc_tree_out)
        else:
            LOGGER.warning("Accessory tree not built (fewer than two samples)")

    # ----------------------------
    # Combined partitions
    #
    # Two samples share a combined partition only if they share BOTH their core
    # partition and their accessory partition, so this is the intersection of the
    # two groupings: every combined partition is a subset of a core partition and
    # of an accessory partition. Samples missing either input (the reference, which
    # is not annotated, or any sample absent from the accessory matrix) are left
    # blank rather than given a singleton, which would imply a real distinction.
    # ----------------------------
    combined_src = {
        sample: (fields['partition'], fields['accessory_partition'])
        for sample, fields in out.items()
        if fields.get('partition') is not None
        and fields.get('accessory_partition') is not None
    }
    if combined_src:
        # sorted so numbering does not depend on dict insertion order
        labels = {key: i + 1 for i, key in enumerate(sorted(set(combined_src.values())))}
        _extend_dict(out, {
            sample: {'combined_partition': labels[key]}
            for sample, key in combined_src.items()
        })
        skipped = len(out) - len(combined_src)
        LOGGER.info(
            "Assigned %d combined partitions across %d samples%s",
            len(labels), len(combined_src),
            f" ({skipped} left blank, missing a core or accessory partition)" if skipped else ""
        )
    else:
        LOGGER.info("Combined partitions not assigned (no sample has both a core and an accessory partition)")

    # ----------------------------
    # Build final summary and write CSV
    # ----------------------------
    summary_out = Path(outdir) / f"{args.prefix}_summary{args.suffix}csv"
    summary_records: List[Dict[str, Any]] = []

    for sample_id, sample_data in out.items():
        rec: Dict[str, Any] = {
            "id": ref_name if sample_id == "Reference" else sample_id,
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
            "accessory_dist_file": {"pane_id": "matrix-3", "file": acc_dist_out},
            "accessory_tree_file": {"pane_id": "tree-2", "file": acc_tree_out},
        }

        # the accessory difference matrix replaces the % Avg Divergence matrix when it is available
        mr_skip = {"mh_dist_file"} if acc_dist_out.exists() else set()
        # panes that are removed completely (pane, settings and file entry) when they are not used
        mr_drop = mr_skip | {"accessory_dist_file", "accessory_tree_file"}

        mr_json = load_json(microreact_template)

        layout = mr_json['panes']['model']['layout']

        for k, v in mr_map.items():
            pane_id, filepath = v.values()
            if k not in mr_skip and _ensure_path(filepath, k):
                _attach_text_file(mr_json, k, filepath)
            else:
                _delete_pane_by_id(layout, pane_id)
                if k in mr_drop:
                    mr_json.get('files', {}).pop(k, None)
                    mr_json.get('matrices', {}).pop(pane_id, None)
                    mr_json.get('trees', {}).pop(pane_id, None)

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
    version = "1.3.0"
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
  - dist_events_all.tsv (optional, pan-dists accessory event distance matrix)
  - *.Rtab (optional, Panaroo gene presence-absence matrix; used only when
    dist_events_all.tsv is absent)
  - template.microreact (optional, for Microreact output)
        """
    )
    parser.add_argument("indir", type=str, help="Input directory containing all required files")
    parser.add_argument("--timestamp", type=str, required=True, help="Run timestamp (epoch)")
    parser.add_argument("--taxa", type=str, required=True, help="Taxonomy/species")
    parser.add_argument("--cluster", type=str, required=True, help="Cluster identifier")
    parser.add_argument("--recomb-masked", action="store_true", help="Recombination masking was performed")
    parser.add_argument("--partition-distance", default=100.0, type=float, help="Distance threshold for core (SNP) tree partitioning, in nucleotides")
    parser.add_argument("--acc-partition-distance", default=2.0, type=float, help="Distance threshold for accessory partitioning, in differing events")
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