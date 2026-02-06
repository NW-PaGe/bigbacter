#!/usr/bin/env python3
"""
bigbacter_report.py
Author: Jared Johnson, jared.johnson@doh.wa.gov

Selects the best reference assembly from multiple FASTA files based on
containment in a global sketch and assembly quality metrics.
"""

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional

import screed
import sourmash as sm


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_contigs(fasta_path: str, min_contig_length: int) -> List[str]:
    """
    Load contigs from FASTA file, filtering by minimum length.
    
    Args:
        fasta_path: Path to FASTA file
        min_contig_length: Minimum contig length to include
        
    Returns:
        List of contig sequences (uppercase)
    """
    if not Path(fasta_path).is_file():
        raise FileNotFoundError(f"File not found: {fasta_path}")
    
    contigs = []
    total_contigs = 0
    
    logger.debug(f"Loading contigs from {fasta_path}")
    
    for record in screed.open(fasta_path):
        total_contigs += 1
        seq = record["sequence"].upper()
        
        if len(seq) >= min_contig_length:
            contigs.append(seq)
    
    filtered_count = total_contigs - len(contigs)
    if filtered_count > 0:
        logger.debug(f"  Filtered {filtered_count}/{total_contigs} contigs below {min_contig_length} bp")
    
    logger.info(f"  Loaded {len(contigs)} contigs ({sum(len(c) for c in contigs):,} bp total)")
    
    return contigs


def build_sketch(contigs: List[str], ksize: int, scaled: int, 
                track_abundance: bool = True) -> sm.MinHash:
    """
    Build MinHash sketch from contig sequences.
    
    Args:
        contigs: List of DNA sequences
        ksize: K-mer size
        scaled: Scaled parameter (downsampling rate)
        track_abundance: Whether to track k-mer abundance
        
    Returns:
        MinHash sketch
    """
    mh = sm.MinHash(
        n=0,
        ksize=ksize,
        scaled=scaled,
        track_abundance=track_abundance,
        seed=42
    )
    
    for seq in contigs:
        mh.add_sequence(seq, force=True)
    
    return mh


def calculate_score(containment: float, length: int, n_contigs: int, 
                   penalty: float) -> float:
    """
    Calculate assembly quality score.
    
    Score = containment × length / (n_contigs ^ penalty)
    
    Higher scores indicate better assemblies (more complete, fewer fragments).
    """
    if n_contigs == 0 or length == 0:
        return float('-inf')
    
    return containment * length / (n_contigs ** penalty)


def get_sample_name(fasta_path: str) -> str:
    """Extract clean sample name from file path."""
    p = Path(fasta_path)
    stem = p.stem
    
    # Handle .fa.gz, .fasta.gz, etc.
    if p.suffix == '.gz':
        stem = Path(stem).stem
    
    return stem


def write_reference_fasta(contigs: List[str], output_path: Path) -> None:
    """Write contigs to gzipped FASTA with simple numeric headers."""
    # Sort contigs by length (smallest to largest)
    sorted_contigs = sorted(contigs, key=len)
    
    logger.info(f"Writing {len(sorted_contigs)} contigs to {output_path}")
    
    with gzip.open(output_path, 'wt', encoding='utf-8') as fh:
        for i, seq in enumerate(sorted_contigs, start=1):
            fh.write(f">{i}\n")
            # Wrap sequence at 60 characters
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j+60] + "\n")


def write_reference_json(name: str, length: int, n_contigs: int, 
                        output_path: Path) -> None:
    """Write reference metadata to JSON file."""
    metadata = {
        "name": name,
        "length": length,
        "n_contigs": n_contigs
    }
    
    logger.info(f"Writing metadata to {output_path}")
    
    with open(output_path, 'w', encoding='utf-8') as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)


def main() -> None:
    """Main entry point."""
    VERSION = "1.0.0"
    
    parser = argparse.ArgumentParser(
        description="Select best reference assembly from multiple FASTA files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        'input',
        nargs='+',
        help='One or more FASTA files to evaluate'
    )
    parser.add_argument(
        '--min-contig-len',
        type=int,
        default=300,
        help='Minimum contig length (bp) to include'
    )
    parser.add_argument(
        '--contig-penalty',
        type=float,
        default=0.5,
        help='Penalize fragmented assemblies (score = d × length / n_contigs^penalty)'
    )
    parser.add_argument(
        '--ksize',
        type=int,
        default=31,
        help='K-mer size for sketching'
    )
    parser.add_argument(
        '--scaled',
        type=int,
        default=100,
        help='Sketch scaled parameter (1/sampling_rate)'
    )
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set log level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Log startup
    logger.info("="*60)
    logger.info(f"BigBacter Report v{VERSION}")
    logger.info(f"Author: Jared Johnson")
    logger.info("="*60)
    logger.info(f"Processing {len(args.input)} assemblies")
    logger.info(f"Parameters: ksize={args.ksize}, scaled={args.scaled}, "
                f"min_contig_len={args.min_contig_len}, penalty={args.contig_penalty}")
    logger.info("")
    
    # Process each assembly
    assemblies = []
    
    for i, fasta_path in enumerate(args.input, 1):
        logger.info(f"[{i}/{len(args.input)}] Processing {fasta_path}")
        
        try:
            # Load and filter contigs
            contigs = load_contigs(fasta_path, args.min_contig_len)
            
            if not contigs:
                logger.warning(f"  No contigs passed filtering - skipping")
                continue
            
            # Build sketch
            logger.debug(f"  Building MinHash sketch")
            mh = build_sketch(contigs, args.ksize, args.scaled)
            
            # Store assembly data
            assemblies.append({
                'name': get_sample_name(fasta_path),
                'path': fasta_path,
                'contigs': contigs,
                'length': sum(len(c) for c in contigs),
                'n_contigs': len(contigs),
                'mh': mh
            })
            
        except Exception as e:
            logger.error(f"  Failed to process {fasta_path}: {e}")
            continue
    
    if not assemblies:
        logger.error("No valid assemblies found - exiting")
        sys.exit(1)
    
    logger.info("")
    logger.info(f"Successfully processed {len(assemblies)} assemblies")
    logger.info("")
    
    # Build global sketch (union of all assemblies)
    logger.info("Building global sketch from all assemblies...")
    global_mh = assemblies[0]['mh'].to_mutable()
    
    for assembly in assemblies[1:]:
        global_mh.merge(assembly['mh'].to_mutable())
    
    logger.info(f"Global sketch contains {len(global_mh)} hashes")
    logger.info("")
    
    # Calculate scores for each assembly
    logger.info("Scoring assemblies:")
    logger.info("-" * 80)
    logger.info(f"{'Sample':<30} {'Length':>12} {'Contigs':>8} {'Contain':>8} {'Score':>12}")
    logger.info("-" * 80)
    
    best_assembly = None
    best_score = float('-inf')
    
    for assembly in assemblies:
        # Calculate containment (what fraction of global diversity is in this assembly)
        containment = global_mh.contained_by(assembly['mh'])
        
        # Calculate score
        score = calculate_score(
            containment,
            assembly['length'],
            assembly['n_contigs'],
            args.contig_penalty
        )
        
        assembly['containment'] = containment
        assembly['score'] = score
        
        # Track best
        if score > best_score:
            best_score = score
            best_assembly = assembly
        
        # Log results
        logger.info(
            f"{assembly['name']:<30} "
            f"{assembly['length']:>12,} "
            f"{assembly['n_contigs']:>8} "
            f"{containment:>8.4f} "
            f"{score:>12.2f}"
        )
    
    logger.info("-" * 80)
    logger.info("")
    
    # Write output files
    logger.info(f"Best assembly: {best_assembly['name']}")
    logger.info(f"  Length: {best_assembly['length']:,} bp")
    logger.info(f"  Contigs: {best_assembly['n_contigs']}")
    logger.info(f"  Containment: {best_assembly['containment']:.4f}")
    logger.info(f"  Score: {best_assembly['score']:.2f}")
    logger.info("")
    
    # Write outputs
    output_fasta = Path("ref.fa.gz")
    output_json = Path("ref.json")
    
    write_reference_fasta(best_assembly['contigs'], output_fasta)
    write_reference_json(
        best_assembly['name'],
        best_assembly['length'],
        best_assembly['n_contigs'],
        output_json
    )
    
    logger.info("")
    logger.info("Output files:")
    logger.info(f"  {output_fasta.resolve()}")
    logger.info(f"  {output_json.resolve()}")
    logger.info("")
    logger.info("Done!")


if __name__ == "__main__":
    main()