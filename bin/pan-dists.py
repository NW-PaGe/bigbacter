#!/usr/bin/env python3
"""
pan-dists.py -- pangenome distances from Panaroo output: accessory event
distances and gene-content distances for every isolate pair.

Model
-----
An EVENT is a cluster of neighbouring accessory genes -- a region of the genome
that is gained or lost as a unit. Events are defined ONCE for the whole
collection, then every isolate pair is scored against that fixed catalogue.

  1. Read gene_presence_absence_roary.csv. A non-empty isolate cell = present.
  2. Classify clusters core (present in >= --core-freq of isolates) or accessory.
  3. Keep only `Genome Fragment` values containing at least one accessory cluster.
  4. Build the event catalogue. Two accessory clusters are linked when they lie
     within N genes of each other on the same contig in at least one isolate
     (--link-mode all requires it in every isolate holding both). Linkage is
     single-linkage, so clusters N+1 apart still join when a third bridges them.
     Connected components of that relation are the events.
  5. Every event is typed by size, and the two types are scored differently:

       MULTI  (>= 2 clusters) -- a block. Compare compositions:
            d == 0   -> identical
            d <= X   -> tolerated, the isolates share this block
            d >  X   -> one event difference
       SINGLE (1 cluster) -- an isolated lone gene with no accessory neighbour
            within N. The tolerance is NEVER applied here, because there is no
            block whose composition could wobble: the gene is either there or
            it is not.
            d == 0   -> identical
            d == 1   -> one event difference

     X therefore only ever acts on multi-cluster events.

  6. Counts are reported three ways, because single-cluster events are the most
     artifact-prone class in draft assemblies and you will want to see the
     number with and without them:
        events_multi  -- block differences only
        events_single -- isolated single-gene differences only
        events_all    -- the sum

Defaults: N = 5, X = 1. Results are written into PANAROO_DIR.

Why positions do not come from `Order within Fragment`
-----------------------------------------------------
Panaroo writes that column as a counter over `for node in
nx.connected_components(G)` -- Python set iteration order, unrelated to genomic
position (generate_output.py). Positions instead come from `clustering_id` in
gene_data.csv, built as <isolate>_<contig>_<gene_index> where gene_index is the
rank of the gene along that contig after sorting by coordinate (prokka.py).

Because linkage pools evidence across isolates, a cluster that is unplaceable in
one isolate (a Panaroo refound gene, `<member>_refound_<n>`, which carries no
positional index) is still placed into its event using the isolates where it was
annotated normally. Only a cluster refound in every isolate holding it stays
unplaced; those form single-cluster events and are reported separately.

Outputs
-------
All files go into one directory (default PANAROO_DIR/pan-dists_N<N>_X<X>), so
no filename prefix is repeated and parameter sweeps never collide.

  events.tsv             the event catalogue: event_type, size, span, members
  pairwise.tsv           one row per pair, every count, split multi/single/all
  detail.tsv             per pair, per non-identical event: d and verdict
  dist_events_all.tsv    distance matrix, multi + single event differences
  dist_events_multi.tsv  distance matrix, multi-cluster event differences only
  dist_genes.tsv         distance matrix, differing gene clusters per pair
  summary.txt            percentiles and caveats
  report.txt             filtering and construction log

The three dist_*.tsv files are square wide-format matrices, ready to compare
against a core SNP matrix.

Usage
-----
  python3 pan-dists.py ./panaroo_out
  python3 pan-dists.py ./panaroo_out -N 8 -X 2
  python3 pan-dists.py ./panaroo_out --always-keep 'bla|mcr|van|tet'
  python3 pan-dists.py ./panaroo_out -o my_results
"""

import argparse
import itertools
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

__version__ = "1.0.0"

ROARY_FIXED_COLUMNS = [
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
    "No. sequences", "Avg sequences per isolate", "Genome Fragment",
    "Order within Fragment", "Accessory Fragment",
    "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc",
]

REFOUND_RE = re.compile(r"_refound_\d+$")
CLUSTERING_ID_RE = re.compile(r"^(\d+)_(\d+)_(\d+)$")
TRAILING_INT_RE = re.compile(r"(\d+)\s*$")


# ----------------------------------------------------------------- input files

def read_roary_csv(path, report):
    df = pd.read_csv(path, low_memory=False, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]

    fixed = [c for c in ROARY_FIXED_COLUMNS if c in df.columns]
    isolate_cols = [c for c in df.columns if c not in fixed]
    if "Gene" not in df.columns:
        sys.exit("No 'Gene' column -- is this a Roary-format CSV?")
    if not isolate_cols:
        sys.exit("No isolate columns found.")

    report.append(f"roary csv: {path}")
    report.append(f"fixed columns found: {len(fixed)}/{len(ROARY_FIXED_COLUMNS)}")
    report.append(f"isolates: {len(isolate_cols)}")
    report.append(f"gene clusters: {len(df)}")
    missing = [c for c in ROARY_FIXED_COLUMNS if c not in df.columns]
    if missing:
        report.append(f"WARNING fixed columns absent: {', '.join(missing)}")
    if "Genome Fragment" not in df.columns:
        report.append("WARNING no 'Genome Fragment' column; fragment filter off")
        df["Genome Fragment"] = "1"

    return df.set_index("Gene", drop=False), isolate_cols


def cell_ids(cell):
    if cell is None:
        return []
    s = str(cell).strip()
    if s == "" or s.lower() == "nan":
        return []
    return [p.strip() for p in s.split(";") if p.strip()]


def read_gene_data(path, report):
    """gene_id -> ((isolate_idx, contig_idx), gene_index), plus scaffold names."""
    if not os.path.exists(path):
        report.append(f"WARNING gene_data.csv not found at {path}")
        return {}, {}, 0

    want = ("gff_file", "scaffold_name", "clustering_id", "annotation_id")
    df = pd.read_csv(path, low_memory=False, dtype=str, keep_default_na=False,
                     usecols=lambda c: c.strip() in want)
    df.columns = [c.strip() for c in df.columns]
    for col in ("clustering_id", "annotation_id", "scaffold_name"):
        if col not in df.columns:
            report.append(f"WARNING gene_data.csv lacks '{col}'")
            return {}, {}, 0

    pos, scaff, n_refound = {}, {}, 0
    for cid, aid, sname in zip(df["clustering_id"], df["annotation_id"],
                               df["scaffold_name"]):
        cid_s, aid_s, sname_s = str(cid).strip(), str(aid).strip(), str(sname).strip()
        m = CLUSTERING_ID_RE.match(cid_s)
        if m:
            entry = ((int(m.group(1)), int(m.group(2))), int(m.group(3)))
            for k in (aid_s, cid_s):
                if k:
                    pos[k] = entry
                    scaff[k] = sname_s
        else:
            if REFOUND_RE.search(cid_s):
                n_refound += 1
            for k in (aid_s, cid_s):
                if k:
                    scaff[k] = sname_s

    report.append(f"gene_data.csv: {len(df)} rows, {len(pos)} ids positioned, "
                  f"{n_refound} refound (no positional id)")
    return pos, scaff, n_refound


def locus_tag_positions(df, isolate_cols, scaff, report):
    pos, n_fail = {}, 0
    for col_i, col in enumerate(isolate_cols):
        for cell in df[col]:
            for gid in cell_ids(cell):
                if gid in pos:
                    continue
                m = TRAILING_INT_RE.search(gid)
                if m and not REFOUND_RE.search(gid):
                    pos[gid] = ((col_i, scaff.get(gid, "")), int(m.group(1)))
                else:
                    n_fail += 1
    report.append(f"locus-tag ordering: {len(pos)} placed, {n_fail} unplaceable")
    if not scaff:
        report.append("WARNING no scaffold info; locus-tag runs may span contig breaks")
    return pos


# -------------------------------------------------------------------- presence

def build_presence(df, isolate_cols):
    pres = pd.DataFrame(False, index=df.index, columns=isolate_cols)
    ids = {col: {} for col in isolate_cols}
    for col in isolate_cols:
        hits = []
        for gene, cell in zip(df.index, df[col]):
            g = cell_ids(cell)
            if g:
                hits.append(gene)
                ids[col][gene] = g
        pres.loc[hits, col] = True
    return pres, ids


def classify(df, pres, args, report):
    n_iso = pres.shape[1]
    counts = pres.sum(axis=1)
    is_core = counts >= args.core_freq * n_iso
    is_acc = (counts >= 1) & (~is_core)
    report.append(f"core clusters (>= {args.core_freq:.2f} of isolates): "
                  f"{int(is_core.sum())}")
    report.append(f"accessory clusters: {int(is_acc.sum())}")

    frag = df["Genome Fragment"].astype(str)
    frags_keep = set(frag[is_acc].unique())
    in_frag = frag.isin(frags_keep)
    report.append(f"genome fragments: {frag.nunique()}, with >=1 accessory: "
                  f"{len(frags_keep)}")
    report.append(f"accessory clusters in retained fragments: "
                  f"{int((is_acc & in_frag).sum())} "
                  f"(dropped {int((is_acc & ~in_frag).sum())})")

    acc = is_acc & in_frag
    if int(acc.sum()) == 0:
        sys.exit("No accessory clusters survive the fragment filter.")
    return acc


# --------------------------------------------------------- the event catalogue

def build_events(acc_clusters, pres, ids, positions, args, report):
    """
    Link accessory clusters that lie within N genes of each other on the same
    contig, pooling evidence across isolates, then take connected components.

    Only consecutive accessory positions need testing: if a and c are within N
    then any b between them is within N of both, so single-linkage via
    connected components recovers the full chain.
    """
    import networkx as nx

    acc_set = set(acc_clusters)
    G = nx.Graph()
    G.add_nodes_from(acc_clusters)

    placed_in = defaultdict(set)      # cluster -> isolates where positioned
    pair_support = defaultdict(int)   # (c1,c2) -> isolates supporting the link
    pair_possible = defaultdict(int)  # (c1,c2) -> isolates holding both, placed

    for iso in pres.columns:
        here = []
        for cluster in acc_clusters:
            if not pres.at[cluster, iso]:
                continue
            hit = None
            for gid in ids[iso].get(cluster, []):
                if gid in positions:
                    hit = positions[gid]
                    break
            if hit is None:
                continue
            (contig, idx) = hit
            here.append((contig, idx, cluster))
            placed_in[cluster].add(iso)

        here.sort(key=lambda e: (str(e[0]), e[1]))
        for (c1, i1, g1), (c2, i2, g2) in zip(here, here[1:]):
            if c1 != c2 or g1 == g2:
                continue
            key = tuple(sorted((g1, g2)))
            pair_possible[key] += 1
            if i2 - i1 <= args.window:
                pair_support[key] += 1

    n_edges = 0
    for key, support in pair_support.items():
        if args.link_mode == "all" and support < pair_possible[key]:
            continue
        G.add_edge(*key)
        n_edges += 1

    unplaced = [c for c in acc_clusters if c not in placed_in]
    report.append(f"linkage edges ({args.link_mode}-isolate, N={args.window}): "
                  f"{n_edges}")
    report.append(f"accessory clusters never positioned in any isolate: "
                  f"{len(unplaced)}")

    comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    events, cluster2event = [], {}
    for k, comp in enumerate(comps):
        eid = f"event_{k:05d}"
        members = sorted(comp)
        for c in members:
            cluster2event[c] = eid
        events.append({"event_id": eid, "n_clusters": len(members),
                       "members": members})
    sizes = [e["n_clusters"] for e in events]
    report.append(f"events defined: {len(events)} "
                  f"(single-cluster {sum(1 for s in sizes if s == 1)}, "
                  f"largest {max(sizes) if sizes else 0})")
    return events, cluster2event, placed_in


def describe_events(events, pres, ids, positions, annot, placed_in):
    """Add a representative contig/span, taken from the isolate holding most members."""
    rows = []
    for ev in events:
        members = ev["members"]
        best = None
        for iso in pres.columns:
            pts = []
            for c in members:
                if not pres.at[c, iso]:
                    continue
                for gid in ids[iso].get(c, []):
                    if gid in positions:
                        pts.append(positions[gid])
                        break
            if not pts:
                continue
            by_contig = defaultdict(list)
            for contig, idx in pts:
                by_contig[contig].append(idx)
            contig, idxs = max(by_contig.items(), key=lambda kv: len(kv[1]))
            if best is None or len(idxs) > best[2]:
                best = (iso, contig, len(idxs), min(idxs), max(idxs))
        n_iso_present = int(pres.loc[members].any(axis=0).sum())
        rows.append({
            "event_id": ev["event_id"],
            "event_type": "single" if ev["n_clusters"] == 1 else "multi",
            "n_clusters": ev["n_clusters"],
            "n_isolates_with_any_member": n_iso_present,
            "ref_isolate": best[0] if best else "",
            "ref_contig": str(best[1]) if best else "",
            "ref_span_start": best[3] if best else "",
            "ref_span_end": best[4] if best else "",
            "ref_span_length": (best[4] - best[3] + 1) if best else "",
            "all_members_positioned": all(c in placed_in for c in members),
            "members": ";".join(members),
            "annotations": ";".join(dict.fromkeys(
                annot.get(m, "") for m in members if annot.get(m, ""))),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------- scoring

def score_pair(a, b, events, pres, annot, always_keep, args):
    """
    Compare compositions event by event, keeping multi-cluster and
    single-cluster events in separate accounts. The tolerance X is applied
    only to multi-cluster events.
    """
    va = pres[a]
    vb = pres[b]
    detail = []
    n_multi = n_single = n_tol = n_forced = 0
    id_multi = id_single = 0
    genes_diff = g_multi = g_single = g_tol = 0

    for ev in events:
        members = ev["members"]
        differing = [c for c in members if bool(va[c]) != bool(vb[c])]
        d = len(differing)
        is_single = ev["n_clusters"] == 1
        etype = "single" if is_single else "multi"
        genes_diff += d

        if d == 0:
            if is_single:
                id_single += 1
            else:
                id_multi += 1
            continue

        if is_single:
            # No block, so no composition to tolerate: X never applies here.
            verdict = "difference_single"
            n_single += 1
            g_single += d
        else:
            forced = bool(always_keep) and any(
                always_keep.search(c) or always_keep.search(annot.get(c, ""))
                for c in differing)
            if d > args.tolerance:
                verdict = "difference_multi"
                n_multi += 1
                g_multi += d
            elif forced:
                verdict = "difference_multi_forced"
                n_multi += 1
                n_forced += 1
                g_multi += d
            else:
                verdict = "tolerated_multi"
                n_tol += 1
                g_tol += d

        in_a = sum(1 for c in differing if bool(va[c]))
        detail.append({
            "isolate_a": a, "isolate_b": b,
            "event_id": ev["event_id"],
            "event_type": etype,
            "n_clusters_in_event": ev["n_clusters"],
            "n_differing": d,
            "frac_differing": round(d / ev["n_clusters"], 4),
            "n_only_in_a": in_a,
            "n_only_in_b": d - in_a,
            "verdict": verdict,
            "differing_clusters": ";".join(differing),
        })

    return detail, {
        "isolate_a": a, "isolate_b": b,
        "events_all": n_multi + n_single,
        "events_multi": n_multi,
        "events_single": n_single,
        "tolerated_multi": n_tol,
        "forced_multi": n_forced,
        "identical_multi": id_multi,
        "identical_single": id_single,
        "genes_different": genes_diff,
        "genes_multi": g_multi,
        "genes_single": g_single,
        "genes_tolerated": g_tol,
    }


# ---------------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(
        prog="pan-dists",
        description="Pangenome distances from Panaroo output: accessory event "
                    "distances and gene-content distances for every pair.")
    p.add_argument("panaroo_dir", metavar="PANAROO_DIR",
                   help="Panaroo output directory; results are written here")
    p.add_argument("-o", "--outdir", default=None,
                   help="output directory "
                        "[PANAROO_DIR/pan-dists_N<N>_X<X>]")
    p.add_argument("--roary-csv", default=None)
    p.add_argument("--gene-data", default=None)
    p.add_argument("-N", "--window", type=int, default=5,
                   help="accessory clusters within N genes on a contig join the "
                        "same event; single-linkage [5]")
    p.add_argument("-X", "--tolerance", type=int, default=1,
                   help="if the compositions of a MULTI-cluster event differ by "
                        "X clusters or fewer, the two isolates are considered "
                        "to share that event. Never applied to single-cluster "
                        "events [1]")
    p.add_argument("--link-mode", default="any", choices=("any", "all"),
                   help="link two clusters if within N in ANY isolate holding "
                        "both (default) or in ALL of them")
    p.add_argument("--core-freq", type=float, default=1.0,
                   help="fraction of isolates at or above which a cluster is "
                        "core [1.0, present in all]")
    p.add_argument("--always-keep", default=None,
                   help="regex on cluster name and annotation; a matching "
                        "difference counts even within tolerance, e.g. "
                        "'bla|mcr|van|tet|erm|qnr'")
    p.add_argument("--order-source", default="clustering-id",
                   choices=("clustering-id", "locus-tag", "fragment-order"))
    p.add_argument("--i-know-this-is-not-an-order", action="store_true",
                   help="required to permit --order-source fragment-order")
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args()

    if not os.path.isdir(args.panaroo_dir):
        sys.exit(f"not a directory: {args.panaroo_dir}")
    roary = args.roary_csv or os.path.join(
        args.panaroo_dir, "gene_presence_absence_roary.csv")
    gene_data = args.gene_data or os.path.join(args.panaroo_dir, "gene_data.csv")
    if not os.path.exists(roary):
        sys.exit(f"missing {roary}")
    if args.outdir is None:
        args.outdir = os.path.join(
            args.panaroo_dir, f"pan-dists_N{args.window}_X{args.tolerance}")
    os.makedirs(args.outdir, exist_ok=True)

    report = [f"N (window) = {args.window}",
              f"X (composition tolerance) = {args.tolerance}",
              f"link mode = {args.link_mode}",
              "tolerance applies to multi-cluster events only",
              f"order source = {args.order_source}", ""]

    df, isolate_cols = read_roary_csv(roary, report)
    pres, ids = build_presence(df, isolate_cols)
    pos_map, scaff, _ = read_gene_data(gene_data, report)

    if args.order_source == "clustering-id":
        positions = pos_map
        if not positions:
            sys.exit("No positions from gene_data.csv. Try --order-source "
                     "locus-tag, or point --gene-data at the right file.")
    elif args.order_source == "locus-tag":
        positions = locus_tag_positions(df, isolate_cols, scaff, report)
    else:
        if not args.i_know_this_is_not_an_order:
            sys.exit(
                "--order-source fragment-order refuses to run.\n"
                "Panaroo writes 'Order within Fragment' as a counter over Python\n"
                "set iteration (generate_output.py), so it is not a genomic\n"
                "order and an N-gene window over it is meaningless. Pass\n"
                "--i-know-this-is-not-an-order to override.")
        report.append("WARNING using 'Order within Fragment' -- NOT a genomic order")
        positions = {}
        oc = pd.to_numeric(df["Order within Fragment"], errors="coerce")
        fc = df["Genome Fragment"].astype(str)
        for col_i, col in enumerate(isolate_cols):
            for gene, cell in zip(df.index, df[col]):
                o = oc.get(gene)
                if pd.isna(o):
                    continue
                for gid in cell_ids(cell):
                    positions.setdefault(gid, ((col_i, fc[gene]), int(o)))

    acc_mask = classify(df, pres, args, report)
    acc_clusters = list(df.index[acc_mask.values])

    annot = (dict(zip(df.index, df["Annotation"].astype(str)))
             if "Annotation" in df.columns else {})
    always_keep = re.compile(args.always_keep, re.I) if args.always_keep else None
    if always_keep:
        report.append(f"--always-keep pattern: {args.always_keep}")

    print(f"defining events over {len(acc_clusters)} accessory clusters "
          f"(N={args.window})")
    events, cluster2event, placed_in = build_events(
        acc_clusters, pres, ids, positions, args, report)

    pairs = list(itertools.combinations(isolate_cols, 2))
    print(f"scoring {len(pairs)} pairs against {len(events)} events "
          f"(X={args.tolerance})")

    detail, rows = [], []
    for a, b in pairs:
        d, summary = score_pair(a, b, events, pres, annot, always_keep, args)
        detail.extend(d)
        rows.append(summary)
    pairwise = pd.DataFrame(rows)

    out = args.outdir

    def path(name):
        return os.path.join(out, name)

    describe_events(events, pres, ids, positions, annot, placed_in).to_csv(
        path("events.tsv"), sep="\t", index=False)
    pairwise.to_csv(path("pairwise.tsv"), sep="\t", index=False)
    DETAIL_COLS = ["isolate_a", "isolate_b", "event_id", "event_type",
                   "n_clusters_in_event", "n_differing", "frac_differing",
                   "n_only_in_a", "n_only_in_b", "verdict", "differing_clusters"]
    (pd.DataFrame(detail) if detail else pd.DataFrame(columns=DETAIL_COLS)
     ).to_csv(path("detail.tsv"), sep="\t", index=False)

    for col, fname in (("events_all", "dist_events_all.tsv"),
                       ("events_multi", "dist_events_multi.tsv"),
                       ("genes_different", "dist_genes.tsv")):
        m = pd.DataFrame(0, index=isolate_cols, columns=isolate_cols, dtype=int)
        for r in pairwise.itertuples(index=False):
            v = getattr(r, col)
            m.at[r.isolate_a, r.isolate_b] = v
            m.at[r.isolate_b, r.isolate_a] = v
        m.to_csv(path(fname), sep="\t")

    n_single_ev = sum(1 for e in events if e["n_clusters"] == 1)
    with open(path("summary.txt"), "w") as fh:
        fh.write(f"N (window) = {args.window}\n"
                 f"X (composition tolerance, multi-cluster events only) = "
                 f"{args.tolerance}\n"
                 f"events in catalogue = {len(events)} "
                 f"({len(events) - n_single_ev} multi, {n_single_ev} single)\n"
                 f"pairs = {len(pairwise)}\n\n")
        for col in ("events_all", "events_multi", "events_single",
                    "tolerated_multi", "genes_different"):
            fh.write(f"{col}\n")
            fh.write(f"  mean {pairwise[col].mean():.3f}  "
                     f"sd {pairwise[col].std():.3f}\n")
            for q in (0, 50, 90, 95, 99, 100):
                fh.write(f"  p{q:<3d} {np.percentile(pairwise[col], q):.3f}\n")
            fh.write("\n")

        tot_all = int(pairwise["events_all"].sum())
        tot_single = int(pairwise["events_single"].sum())
        frac = (tot_single / tot_all) if tot_all else 0.0
        fh.write(f"single-cluster events as a share of events_all: "
                 f"{tot_single}/{tot_all} = {frac:.1%}\n")
        if frac > 0.25:
            fh.write("  That is a large share. Isolated single-gene calls are\n"
                     "  the most artifact-prone class in draft assemblies\n"
                     "  (contig-boundary truncation, missed start codons,\n"
                     "  frameshifts splitting one gene into two clusters).\n"
                     "  Prefer events_multi as the headline figure and treat\n"
                     "  events_single as a separate, noisier signal.\n")
        n_tol = int(pairwise["tolerated_multi"].sum())
        fh.write(f"\nmulti-cluster comparisons tolerated across all pairs: "
                 f"{n_tol}\n")
        if n_tol:
            fh.write("  each is a row in detail.tsv with verdict\n"
                     "  'tolerated_multi' -- read them before trusting counts.\n")
        if len(pairwise) < 30:
            fh.write(f"\nNOTE only {len(pairwise)} pairs. Percentiles are not a\n"
                     "distribution at this sample size -- read detail.tsv\n"
                     "and judge each event rather than applying a cutoff.\n")

    with open(path("report.txt"), "w") as fh:
        fh.write("\n".join(report) + "\n")

    print("\n".join(report))
    print(f"\nwrote to {out}/ :")
    for name in ("events.tsv", "pairwise.tsv", "detail.tsv",
                 "dist_events_all.tsv", "dist_events_multi.tsv",
                 "dist_genes.tsv", "summary.txt", "report.txt"):
        print(f"      {name}")


if __name__ == "__main__":
    main()
