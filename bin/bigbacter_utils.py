
import logging
import os
import sys
import screed
import re
from pathlib import Path

def get_assembly_stem(path):
    path = Path(path)
    stem = path.name
    suffixes = path.suffixes

    if not suffixes:
        return stem

    if suffixes[-1] == ".gz":
        drop_list = suffixes[-2:]   # drop ".fastq" and ".gz"
    else:
        drop_list = suffixes[-1:]   # drop last suffix only

    for s in drop_list:
        stem = stem.removesuffix(s)

    return stem

def safe_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)

def get_ref_name(rec, cache):
    def assign_name(rec):
        # option 1 - use existing name
        name = rec.get('name') or rec.get('metadata', {}).get('name')
        if name:
            return str(name)

        # option 2 - name reference using taxon, segment, variant fields
        taxon   = rec.get('taxon') or rec.get('metadata', {}).get('taxon')
        segment = rec.get('segment') or rec.get('metadata', {}).get('segment')
        variant = rec.get('variant') or rec.get('metadata', {}).get('variant')

        if taxon and segment and variant:
            return safe_filename(f"{str(taxon)}-{str(segment)}-{str(variant)}")

        # option 3 - name reference using basename field
        basename = rec.get('basename')
        if basename:
            return basename

        # option 4 - provide generic name
        return "Reference"

    name = safe_filename(assign_name(rec))

    if name in cache:
        cache[name] += 1
        name = f"{name}_{str(cache[name])}"
    else:
        cache[name] = 0

    return name

def logging_config(log_level='INFO'):
    # Get script name
    script_name = os.path.basename(sys.argv[0]).replace('.py', '')
    LOGGER = logging.getLogger(script_name)

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    return LOGGER


def load_single_fasta_record(path: str, contig: str | None = None):
    """
    Load a single FASTA/FASTQ record.

    - If `contig` is provided: return only that record; error if not found.
    - If `contig` is None or "", return the first record and ignore the rest.
    """
    rec = None
    found = False

    for r in screed.open(path):
        name = r.name.split()[0]  # normalize (strip after whitespace)

        if contig:  # contig name supplied → select only that one
            if name == contig:
                return r  # return EXACT record immediately
        else:  # no contig supplied → take the first and ignore the rest
            return r

    # After loop → check failure conditions
    if contig:
        raise ValueError(f"Contig '{contig}' not found in file: {path}")
    else:
        raise ValueError(f"File contains no records: {path}")
