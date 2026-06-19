"""
Download and prepare GSE104406 for a young-versus-old differential expression
analysis of human hematopoietic stem cells (Adelman et al. 2019, Cancer Discov).

Design: 10 young donors (ages 18 to 30) vs 10 aged donors (ages 65 to 75),
FACS-sorted bone-marrow HSCs, bulk RNA-seq. A clean balanced two-group design.

We use NCBI's uniformly processed RNA-seq counts (raw gene-level counts produced
by the GEO/SRA pipeline), so no read alignment is needed. The sample-to-group
labels come from the GEO series matrix, parsed at runtime rather than hardcoded.

Sources:
  raw counts:  https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html
  series GSE104406: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104406

Outputs:
  data/counts.csv     genes (Entrez GeneID) x samples, integer raw counts
  data/meta.csv       sample, condition (young or old)
  data/gene_map.csv   gene_id, symbol
"""

import gzip
import io
import os
import re
import sys
import urllib.request

import pandas as pd

ACC = "GSE104406"
COUNTS_URL = ("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
              "&acc=GSE104406&format=file&file=GSE104406_raw_counts_GRCh38.p13_NCBI.tsv.gz")
ANNOT_URL = ("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
             "&acc=GSE104406&format=file&file=Human.GRCh38.p13.annot.tsv.gz")
SERIES_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/"
              "GSE104406/matrix/GSE104406_series_matrix.txt.gz")

# Drop genes with fewer than this many reads summed across all samples. Standard
# low-count pre-filter; it removes noise and speeds up the model.
MIN_TOTAL_COUNT = 10

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")


def fetch(url, dest):
    """Download url to dest unless it already exists."""
    if os.path.exists(dest):
        return dest
    print("[download] fetching {}".format(os.path.basename(dest)))
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def classify(text):
    """Assign a sample to 'young' or 'old' from its title and characteristics."""
    t = text.lower()
    if "young" in t:
        return "young"
    if "old" in t or "aged" in t or "elderly" in t:
        return "old"
    # Fallback: read an explicit age number. Young 18 to 30, aged 65 to 75, so a
    # 40 / 60 split is safe and never lands between the two groups.
    m = re.search(r"age[^0-9]*([0-9]{1,3})", t)
    if m:
        age = int(m.group(1))
        if age <= 40:
            return "young"
        if age >= 60:
            return "old"
    return None


def parse_series_groups(series_path):
    """Parse the GEO series matrix and return {GSM: 'young'/'old'}."""
    gsms, titles = [], []
    characteristics = []  # list of per-line lists, each aligned to gsms
    with gzip.open(series_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsms = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics"):
                characteristics.append(
                    [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]])
    if not gsms:
        sys.exit("Could not parse sample accessions from the series matrix.")

    groups = {}
    for i, gsm in enumerate(gsms):
        parts = []
        if i < len(titles):
            parts.append(titles[i])
        for row in characteristics:
            if i < len(row):
                parts.append(row[i])
        groups[gsm] = classify(" ".join(parts))
    return groups


def main():
    counts_path = fetch(COUNTS_URL, os.path.join(RAW_DIR, ACC + "_raw_counts.tsv.gz"))
    annot_path = fetch(ANNOT_URL, os.path.join(RAW_DIR, "Human.GRCh38.p13.annot.tsv.gz"))
    series_path = fetch(SERIES_URL, os.path.join(RAW_DIR, ACC + "_series_matrix.txt.gz"))

    # Raw counts: first column is GeneID, remaining columns are GSM samples.
    print("[download] reading counts matrix")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index.name = "gene_id"
    print("[download] counts: {} genes x {} samples".format(*counts.shape))

    # Map each sample (GSM) to young or old from the series metadata.
    groups = parse_series_groups(series_path)
    sample_group = {s: groups.get(s) for s in counts.columns}
    unlabeled = [s for s, g in sample_group.items() if g is None]
    if unlabeled:
        print("[download] WARNING: could not label these samples:", unlabeled)
    meta = pd.DataFrame(
        {"sample": list(sample_group.keys()), "condition": list(sample_group.values())}
    ).dropna(subset=["condition"]).set_index("sample")

    counts = counts[meta.index]  # keep only labeled samples, aligned order
    n_young = int((meta["condition"] == "young").sum())
    n_old = int((meta["condition"] == "old").sum())
    print("[download] labeled samples: {} young, {} old".format(n_young, n_old))
    if n_young < 2 or n_old < 2:
        sys.exit("Need at least 2 samples per group; check the parsed labels above.")

    # Low-count pre-filter: drop genes with too few reads across all samples.
    keep = counts.sum(axis=1) >= MIN_TOTAL_COUNT
    counts = counts[keep]
    print("[download] kept {} genes after low-count filter (>= {} total reads)".format(
        counts.shape[0], MIN_TOTAL_COUNT))

    # Gene annotation: GeneID -> Symbol, for readable result tables and plots.
    print("[download] reading gene annotation")
    annot = pd.read_csv(annot_path, sep="\t", index_col=0)
    sym_col = "Symbol" if "Symbol" in annot.columns else annot.columns[0]
    gene_map = annot.loc[annot.index.intersection(counts.index), [sym_col]].copy()
    gene_map.columns = ["symbol"]
    gene_map.index.name = "gene_id"

    os.makedirs(DATA_DIR, exist_ok=True)
    counts.to_csv(os.path.join(DATA_DIR, "counts.csv"))
    meta.to_csv(os.path.join(DATA_DIR, "meta.csv"))
    gene_map.to_csv(os.path.join(DATA_DIR, "gene_map.csv"))
    print("[download] wrote data/counts.csv, data/meta.csv, data/gene_map.csv")


if __name__ == "__main__":
    main()
