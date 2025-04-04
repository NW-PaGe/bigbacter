#!/usr/bin/env python3

# bb-cluster.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version = 1.0

# Bioinformatics modules
import sourmash
from sourmash import SourmashSignature, save_signatures, load_one_signature
import screed

# Scientific computing and data processing
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS

# General utilities
import random
import csv
import datetime
import argparse
import itertools
import re
import os
import shutil

# Visualization
import matplotlib.pyplot as plt

# ----- ARGS -----#
# Initialize the parser
parser = argparse.ArgumentParser(description="Genomic clustering tool for population-level surviellance.")
# Add arguments
parser.add_argument(
    "--query", dest="query", nargs='+', required=True,
    help="Path(s) to query FASTA(s)."
)
parser.add_argument(
    "--sigs", dest="sigs", nargs='+', required=False,
    help="Path(s) to existing signatures."
)
parser.add_argument(
    "--clusters", dest="clusters", required=False,
    help="Path to file containing cluster information for signatures."
)
parser.add_argument(
    "--dist", dest="dist", default=0.1, type=float,
    help="Distance thresholds to create clusters. Default is 0.1."
)
parser.add_argument(
    "--ksize", dest="ksize", default=31, type=int,
    help="k-mer size used by Sourmash. Default is 31."
)
parser.add_argument(
    "--outdir", dest="outdir", default='./', type=str,
    help="Path to output directory. Default is current directory."
)
parser.add_argument(
    "--summarize", dest="summarize", action='store_true',
    help="Create a summary between all isolates in a database. Can be slow for large databases."
)
parser.add_argument(
    "--no-update", dest="no_update", action='store_true',
    help="Do not include signatures supplied via '--sigs' in the output."
)
parser.add_argument(
    "--version", action='version', version=f'{version}',
    help="Show program version and exit."
)
# Parse the arguments
args = parser.parse_args()

#----- START MESSAGE -----#
start = f"""
bb-cluster.py v{version}

Written by Jared Johnson
jared.johnson@doh.wa.gov

"""
print(start)

#----- FUNCTIONS -----#
def table2term(data, header):
    """
    Prints a list of lists to the terminal with equal padding and header

    Parameters:
    - data (list of lists): The input list of sublists.
    - header (str): Header that will be printed above the data

    Returns:
    - Prints result to terminal
    """
    padded_list = []
    
    # Determine the length of the longest element across all sublists
    max_col_len = max([len(str(item)) for sublist in data for item in sublist]) + 2
    for sublist in data:
        # Pad each element in the sublist to match the maximum length
        padded_line = f'{" "*4}'.join([ str(i).ljust(max_col_len) for i in sublist])
        padded_list.append(padded_line)
    max_row_len = max([len(str(row)) for row in padded_list])
    padded_result = '\n'.join(padded_list)
    if len(str(header)) > max_row_len:
        header_pad = 0
    else:
        header_pad = int(round((max_row_len - len(header)) / 2, 0))
    output = f"\n{'-' * max_row_len}\n{' ' * header_pad}{header}\n{'-' * max_row_len}\n{padded_result}\n{'-' * max_row_len}\n"
    print(output, flush=True)

def dictToCsv(dicts, filename, cols):
    """
    Writes a list of dictionaries to a CSV file.

    Parameters:
    - dicts (list of dict): Data to write to CSV.
    - filename (str): The name of the output file.
    - cols (list of str): Column headers for the CSV file.
    """
    # Open the file for writing
    with open(filename, 'w', newline='') as file:
        # Create a CSV writer with specified column headers
        writer = csv.DictWriter(file, fieldnames=cols, quoting=csv.QUOTE_ALL)

        # Write the header row
        writer.writeheader()

        # Write the data rows
        writer.writerows(dicts)

    print(f"Data written to {filename}")


def computeDistance(pair, data):
    """
    Compute the ANI (Average Nucleotide Identity) distance between a pair of items.

    Parameters:
    - pair (tuple): A tuple containing the IDs of two items (id1, id2).
    - data (dict): A dictionary where each ID maps to metadata, 
                   including a 'mh' object with an avg_containment method.

    Returns:
    - dict: A dictionary with the IDs of the items and the computed distance.
    """
    # Unpack the pair of IDs
    id1, id2 = pair

    # Calculate the average containment value
    avg_containment = data[id1]['mh'].avg_containment(data[id2]['mh'])

    # Return the computed distance
    return {
        'id1': id1,
        'id2': id2,
        'dist': 1 - avg_containment
    }


def computeMatrix(data):
    """
    Computes matrix based on pairwise ANI distances.

    Parameters:
    - data (dict): A dictionary where each key maps to metadata required for distance computation.

    Returns:
    - mat (numpy.ndarray): A square matrix of distances.
    - ids (list): A sorted list of unique IDs.
    - dist_long (list of dict): Pairwise comparisons in long format.
    """
    dist_long = []

    # Get unique IDs and sort them
    ids = sorted(set(data.keys()))

    # Generate all pairwise comparisons
    pairs = list(itertools.combinations(ids, 2))

    # Compute distances for each pair
    for pair in pairs:
        dist_long.append(computeDistance(pair, data))

    # Write pairwise comparisons to CSV file
    dictToCsv(
        dist_long, 
        os.path.join(args.outdir, 'dist.new.csv'), 
        ["id1", "id2", "dist"]
    )

    # Extract unique IDs involved in comparisons
    ids = sorted(set(d['id1'] for d in dist_long).union(d['id2'] for d in dist_long))

    # Initialize the distance matrix with zeros
    mat = [[0] * len(ids) for _ in range(len(ids))]

    # Populate the matrix with distances
    for entry in dist_long:
        i = ids.index(entry['id1'])
        j = ids.index(entry['id2'])
        mat[i][j] = entry['dist']
        mat[j][i] = entry['dist']

    # Convert the matrix to a numpy array and finalize formatting
    mat = np.array(mat)
    np.fill_diagonal(mat, 0)  # Ensure diagonal is zero
    mat = squareform(mat)     # Convert to condensed distance matrix format

    return mat, ids, dist_long


def createClusters(data, threshold, start, stage):
    """
    Assigns clusters to items using a distance threshold.

    Parameters:
    - data (dict): A dictionary containing item metadata.
    - threshold (float): Distance threshold for clustering.
    - start (int): Starting index for cluster numbering.
    - stage (str): The stage name for plotting and saving dendrograms.

    Returns:
    - result (dict): A dictionary mapping items to cluster IDs.
    - max_cluster (int): The maximum cluster index.
    - dist_long (list of dict): Pairwise distance information.
    """
    # Handle edge case: single-item data
    if len(data.items()) == 1:
        max_cluster = int(start) + 1
        return {next(iter(data)): max_cluster}, max_cluster, []

    result = {}

    # Step 1: Compute pairwise distances
    mat, ids, dist_long = computeMatrix(data)

    # Step 2: Perform hierarchical clustering
    Z = linkage(mat, method='complete')

    # Step 3: Assign clusters based on the distance threshold
    clusters = cut_tree(Z, height=threshold).flatten()
    for index, item_id in enumerate(ids):
        result[item_id] = int(clusters[index])

    # Step 4: Plot and save dendrogram
    try:
        plt.figure(figsize=(10, 5))
        dendro = dendrogram(Z, labels=ids, orientation='left')
        plt.axvline(x=threshold, color='r', linestyle='--')
        plt.title(stage)
        plt.xlabel('Samples')
        plt.ylabel('Distance')
        plt.savefig(os.path.join(args.outdir, f"{stage}.jpg"), format='png')
        plt.close()
    except Exception as e:
        print(f"Dendrogram not made: {e}", flush=True)

    # Step 5: Merge overlapping clusters
    clusts2merge = {}
    for pair in dist_long:
        id1, id2, dist = pair['id1'], pair['id2'], pair['dist']
        clust1, clust2 = result[id1], result[id2]
        if dist <= threshold and clust1 != clust2:
            clusts2merge.setdefault(clust1, []).append(clust2)
            clusts2merge.setdefault(clust2, []).append(clust1)

    # Step 6: Resolve indirect overlaps
    merge_tbl = {}
    for key, value in clusts2merge.items():
        if key not in merge_tbl:
            overlaps = [key]
            new_overlaps = list(set(value))
            while new_overlaps:
                overlap_update = []
                for cluster in new_overlaps:
                    overlap_update.extend(c for c in clusts2merge.get(cluster, []) if c not in overlaps)
                    overlaps.append(cluster)
                new_overlaps = list(set(overlap_update))
            for cluster in overlaps:
                merge_tbl[cluster] = min(overlaps)

    # Step 7: Update and rename clusters
    final_clusts = []
    for item_id, cluster_id in result.items():
        result[item_id] = merge_tbl.get(cluster_id, cluster_id)
        final_clusts.append(result[item_id])
    final_clusts = sorted(set(final_clusts))

    max_cluster = len(final_clusts) + int(start)
    clust_tbl = {clust: index + int(start) + 1 for index, clust in enumerate(final_clusts)}
    for item_id, cluster_id in result.items():
        result[item_id] = clust_tbl[cluster_id]

    return result, max_cluster, dist_long


def selectReps(minhash, clusters):
    """
    Select representative sequences for each cluster based on minimal average intra-cluster distance.

    Parameters:
    - minhash (dict): A dictionary where keys are isolate IDs, and values are metadata, including 'mh' objects.
    - clusters (dict): A dictionary mapping isolates to their cluster IDs.

    Returns:
    - reps (dict): A dictionary containing the representative sequence and metadata for each cluster.
    """
    data = {}
    reps = {}

    # Organize clusters and associate isolates with their minhash objects
    for isolate, cluster_id in clusters.items():
        if cluster_id not in data:
            data[cluster_id] = []
            reps[cluster_id] = {}
        data[cluster_id].append([isolate, minhash[isolate]['mh']])

    # Compute representatives for each cluster
    for cluster_id, isolates in data.items():
        if len(isolates) == 1:
            reps[cluster_id]['min_avg_dist'] = 0
            reps[cluster_id]['rep'] = isolates[0][0]
            reps[cluster_id]['mh'] = minhash[isolates[0][0]]['mh']
            continue

        for isolate_1 in isolates:
            dists = [
                1 - isolate_1[1].avg_containment(isolate_2[1])
                for isolate_2 in isolates if isolate_1[0] != isolate_2[0]
            ]
            avg_dist = sum(dists) / len(dists)

            if 'min_avg_dist' not in reps[cluster_id] or avg_dist < reps[cluster_id]['min_avg_dist']:
                reps[cluster_id]['min_avg_dist'] = avg_dist
                reps[cluster_id]['rep'] = isolate_1[0]
                reps[cluster_id]['mh'] = minhash[isolate_1[0]]['mh']

    # Print representative sequences to the console
    summary_lines = [['Cluster', 'Isolate', 'Avg Intra Dist']]
    for cluster_id, cluster_data in reps.items():
        summary_lines.append([cluster_id, cluster_data['rep'], round(cluster_data['min_avg_dist'], 3)])
    table2term(
        summary_lines,
        "CLUSTER REPRESENTATIVES"
    )

    return reps

def assignClusters(reps, new, threshold):
    """
    Assign clusters to new items based on existing cluster representatives.

    Parameters:
    - reps (dict): A dictionary of existing cluster representatives. 
                   Keys are cluster IDs, values are metadata including 'mh' objects.
    - new (dict): A dictionary of new items to be clustered. 
                  Keys are item IDs, values are metadata including 'mh' objects.
    - threshold (float): The distance threshold for assigning items to clusters.

    Returns:
    - clusters (dict): A dictionary mapping new item IDs to their assigned cluster IDs.
    - remainder (list): A list of item IDs that could not be assigned to any cluster.
    """
    remainder = []  # Items that couldn't be assigned to a cluster
    clusters = {}   # Mapping of item IDs to cluster IDs

    # Iterate through new items
    for new_id, new_metadata in new.items():
        assigned_cluster = {}

        # Compare each new item to the existing cluster representatives
        for cluster_id, rep_metadata in reps.items():
            dist = 1 - new_metadata['mh'].avg_containment(rep_metadata['mh'])
            if dist < threshold:
                assigned_cluster[cluster_id] = dist

        # Assign the item to the closest cluster if within the threshold
        if assigned_cluster:
            clusters[new_id] = min(assigned_cluster, key=assigned_cluster.get)
        else:
            # Otherwise, add to the remainder
            remainder.append(new_id)

    return clusters, remainder


def summarize(data, clusters):
    """
    Creates a full database summary by computing cluster statistics, pairwise distances, and visualizing with MDS.

    Parameters:
    - data (dict): A dictionary containing metadata for each sample.
    - clusters (dict): A dictionary mapping sample IDs to cluster IDs.

    Outputs:
    - Saves CSV files containing summary data and pairwise comparisons.
    - Saves MDS projection as an image file.
    - Prints cluster statistics to the console.
    """
    dist_long = []
    clustCounts = {}

    # Step 1: Count samples per cluster
    for sample_id, cluster_id in clusters.items():
        clustCounts[str(cluster_id)] = clustCounts.get(str(cluster_id), 0) + 1

    # Step 2: Compute pairwise distances
    ids = sorted(data.keys())
    pairs = itertools.combinations(ids, 2)
    for pair in pairs:
        dist_long.append({
            **computeDistance(pair, data),
            **{'cluster1': clusters[pair[0]], 'cluster2': clusters[pair[1]]}
        })

    dictToCsv(
        dist_long, 
        os.path.join(args.outdir, 'dist.all.csv'), 
        ["id1", "id2", "dist", "cluster1", "cluster2"]
    )

    # Step 3: Construct distance matrix
    ids = sorted(set(d['id1'] for d in dist_long).union(d['id2'] for d in dist_long))
    mat = np.zeros((len(ids), len(ids)))
    for entry in dist_long:
        i, j = ids.index(entry['id1']), ids.index(entry['id2'])
        mat[i][j] = mat[j][i] = entry['dist']
    np.fill_diagonal(mat, 0)

    # Step 4: Summarize cluster relationships
    clust2clust = {}
    for pair in dist_long:
        id1, id2 = str(pair['id1']), str(pair['id2'])
        c1, c2, dist = str(pair['cluster1']), str(pair['cluster2']), pair['dist']
        if id1 != id2:
            clust2clust.setdefault(c1, {}).setdefault(c2, []).append(dist)
            clust2clust.setdefault(c2, {}).setdefault(c1, []).append(dist)

    clusts = set(str(clusters[id]) for id in ids)
    fullSmry = []
    summary_lines = [["Cluster", "n", "Neighbor", "Intra Min", "Intra Max", "Inter Min", "Inter Max"]]
    for c1 in clusts:
        subsummary = [c1, clustCounts[c1], '', 0, 0, 0, 0]
        for c2 in clusts:
            dists = clust2clust.get(c1, {}).get(c2, [0])
            fullSmry.append({'cluster1': c1, 'cluster2': c2, 'min': min(dists), 'max': max(dists), 'mean': sum(dists) / len(dists)})
            if c1 == c2:
                subsummary[3], subsummary[4] = round(min(dists), 3), round(max(dists), 3)
            else:
                d_min, d_max = min(dists), max(dists)
                if d_min < subsummary[5] or subsummary[2] == '':
                    subsummary[5], subsummary[2] = round(d_min, 3), c2
                elif d_min == subsummary[5]:
                    subsummary[2] += f"; {c2}"
                if d_max > subsummary[6]:
                    subsummary[6] = round(d_max, 3)
        summary_lines.append(list(map(str, subsummary)))

    dictToCsv(
        fullSmry, 
        os.path.join(args.outdir, 'summary.full.csv'), 
        ["cluster1", "cluster2", "min", "max", "mean"]
    )

    summary_string = "\n".join([ ','.join(line) for line in summary_lines ])
    with open(os.path.join(args.outdir, 'summary.simple.csv'), "w") as file:
        file.write(summary_string)

    table2term(
        summary_lines,
        "RUN SUMMARY"
    )

    # Step 5: Visualize clusters using MDS
    try:
        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, normalized_stress="auto")
        embedded_coordinates = mds.fit_transform(mat)
        x, y = embedded_coordinates[:, 0], embedded_coordinates[:, 1]

        labels = [str(c) for c in clusts]
        for cluster, label in zip(clusts, labels):
            indices = [i for i, clust in enumerate(clusts) if clust == cluster]
            plt.scatter([x[i] for i in indices], [y[i] for i in indices], label=label)
        
        plt.title("MDS Projection")
        plt.xlabel("Dimension 1")
        plt.ylabel("Dimension 2")
        plt.legend()
        plt.savefig(os.path.join(args.outdir, 'MDS.jpg'), dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"MDS not made: {e}", flush=True)

#----- MAIN -----#
# Configuration and setup
## Create output directory, if needed
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

## Validate input: Signatures and clusters must be supplied together
if (args.clusters and not args.sigs) or (args.sigs and not args.clusters):
    sys.exit("ERROR: Cluster information and signatures must be supplied together.")

# Create signature directory
sig_dir = os.path.join(args.outdir, 'sig/')
if not os.path.exists(sig_dir):
    os.makedirs(sig_dir)

# Generate MinHash signatures for query genomes
mh_query = {}
for genome in args.query:
    name = os.path.splitext(os.path.basename(genome))[0]
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1000)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=name)
    with open(os.path.join(sig_dir, f'{name}.sig'), 'wt') as fp:
        save_signatures([sig], fp)
    mh_query[name] = {'mh': mh}

print(f'{datetime.datetime.now()}: Sketched {len(mh_query)} query sequences.', flush=True)

# Handle existing clusters and signatures
if args.clusters and args.sigs:
    # Load cluster data
    old_clusters = {}
    with open(args.clusters, mode='r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            old_clusters[row[0]] = row[1]

    # Load existing signatures
    mh_sig = {}
    for sig in args.sigs:
        sig_name = os.path.splitext(os.path.basename(sig))[0]
        if sig_name in old_clusters:
            mh = load_one_signature(sig)
            mh_sig[sig_name] = {'mh': mh.minhash}
            if not args.no_update:
                shutil.copy(sig, os.path.join(sig_dir, os.path.basename(sig)))
        else:
            print(f"{sig_name} not found in supplied clusters - it will be skipped!", flush=True)

    # Verify at least one valid signature exists
    if not mh_sig:
        sys.exit("ERROR: Signature sample names do not match cluster sample names!")

    # Assign new isolates to existing clusters
    print(f'{datetime.datetime.now()}: Assigning new isolates to existing clusters.', flush=True)
    mh_reps = selectReps(mh_sig, old_clusters)
    assigned_clusters, remainder = assignClusters(mh_reps, mh_query, args.dist)
    print(f'{datetime.datetime.now()}: {len(assigned_clusters)} isolates assigned.', flush=True)

    # Handle remaining unassigned isolates
    mh_rem = {s: mh_query[s] for s in remainder}
    max_cluster = max(mh_reps.keys()) if mh_rem else 0
else:
    mh_rem = mh_query
    mh_sig = {}
    old_clusters = {}
    assigned_clusters = {}
    max_cluster = 0

# Create new clusters for unassigned isolates
if mh_rem:
    print(f'{datetime.datetime.now()}: Creating new clusters.', flush=True)
    new_clusters, new_max_cluster, dist_long = createClusters(mh_rem, args.dist, max_cluster, 'new_clusters')
    print(f'{datetime.datetime.now()}: {int(new_max_cluster) - int(max_cluster)} clusters created.', flush=True)
else:
    new_clusters = {}
    new_max_cluster = max_cluster

# Combine all clusters
clusters = {**old_clusters, **assigned_clusters, **new_clusters}

# Save cluster assignments to a file
with open(os.path.join(args.outdir, 'clusters.csv'), mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["name", "cluster"])
    for key, value in clusters.items():
        writer.writerow([key, value])

# Create run summary, if requested
if args.summarize:
    print(f'{datetime.datetime.now()}: Creating run summary.', flush=True)
    summarize({**mh_query, **mh_sig}, clusters)

# Print completion message
print(f"\nResults saved to {args.outdir}", flush=True)
