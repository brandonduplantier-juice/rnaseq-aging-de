"""
Download and prepare a young-versus-old human RNA-seq dataset from GEO for
differential expression.

Target: GSE104704 (Nativio et al.), human prefrontal cortex, healthy young vs
healthy aged brains (a third diseased/Alzheimer group is present and dropped, so
this is a clean healthy-aging comparison).

GEO is inconsistent about how counts are stored, so this loader tries two sources
in order and reports what it finds:
  1. NCBI's uniformly processed raw counts (fast path, columns are GSM IDs).
  2. The series' own supplementary directory, auto-discovering a raw-counts file.
If neither yields a usable counts matrix, it prints the available supplementary
files and exits so we can point it at the right file.

Sample group labels (young vs old) are parsed from the GEO series matrix at
runtime and printed for inspection, never hardcoded.

Outputs:
  data/counts.csv     genes x samples, integer raw counts
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

ACC = "GSE104704"
MIN_TOTAL_COUNT = 10        # drop genes with fewer than this many reads in total

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")


def stub(acc):
    """GSE104704 -> GSE104nnn (GEO FTP directory bucket)."""
    digits = acc[3:]
    return "GSE" + digits[:-3] + "nnn"


NCBI_COUNTS_URL = ("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
                   "&acc={a}&format=file&file={a}_raw_counts_GRCh38.p13_NCBI.tsv.gz").format(a=ACC)
NCBI_ANNOT_URL = ("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
                  "&acc={a}&format=file&file=Human.GRCh38.p13.annot.tsv.gz").format(a=ACC)
SUPPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{s}/{a}/suppl/".format(s=stub(ACC), a=ACC)
SERIES_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/{s}/{a}/matrix/"
              "{a}_series_matrix.txt.gz").format(s=stub(ACC), a=ACC)


def http_get(url):
    """Return raw bytes for a URL, or None on HTTP error."""
    try:
        with urllib.request.urlopen(url) as r:
            return r.read()
    except urllib.error.HTTPError as e:
        print("[download]   HTTP {} for {}".format(e.code, url))
        return None


def save(data, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(data)
    return path


def classify(text):
    """Map a sample to young, old, or None (None = drop, e.g. diseased)."""
    t = text.lower()
    # This dataset includes Alzheimer's brains; we compare healthy young vs old,
    # so exclude any diseased sample first.
    if "alzheimer" in t or "diseas" in t or re.search(r"\bad\b", t):
        return None
    if "young" in t:
        return "young"
    if "old" in t or "aged" in t or "elderly" in t:
        return "old"
    m = re.search(r"age[^0-9]*([0-9]{1,3})", t)
    if m:
        age = int(m.group(1))
        if age <= 45:
            return "young"
        if age >= 60:
            return "old"
    return None


def parse_series(series_bytes):
    """Return (gsm_to_group, gsm_to_title) parsed from the series matrix bytes."""
    gsms, titles, chars = [], [], []
    with gzip.open(io.BytesIO(series_bytes), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            if line.startswith("!Sample_geo_accession"):
                gsms = cells[1:]
            elif line.startswith("!Sample_title"):
                titles = cells[1:]
            elif line.startswith("!Sample_characteristics"):
                chars.append(cells[1:])
    gsm_to_group, gsm_to_title = {}, {}
    print("[download] parsed {} samples from series matrix:".format(len(gsms)))
    for i, gsm in enumerate(gsms):
        title = titles[i] if i < len(titles) else ""
        text_parts = [title] + [row[i] for row in chars if i < len(row)]
        grp = classify(" ".join(text_parts))
        gsm_to_group[gsm] = grp
        gsm_to_title[gsm] = title
        print("  {}  {:<45}  -> {}".format(gsm, title[:45], grp))
    return gsm_to_group, gsm_to_title


def find_counts_file(listing_html):
    """Pick a raw-counts file name from a suppl directory HTML listing."""
    names = re.findall(r'href="([^":?/][^"]*\.(?:gz|txt|tsv|csv))"', listing_html)
    names = sorted(set(names))
    skip = ("norm", "fpkm", "tpm", "rpkm", "cpm")
    for pattern in (r"raw[_.-]?count", r"\bcount", r"gene[_.-]?count"):
        for n in names:
            nl = n.lower()
            if re.search(pattern, nl) and not any(k in nl for k in skip):
                return n, names
    return None, names


def read_counts_bytes(data, fname):
    """Read a gzip or plain counts file (bytes) into a DataFrame, gene rows."""
    raw = gzip.decompress(data) if fname.endswith(".gz") else data
    text = raw.decode("utf-8", "replace")
    sep = "\t" if (text[:5000].count("\t") >= text[:5000].count(",")) else ","
    df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
    df.index.name = "gene_id"
    return df


def attach_symbols(counts):
    """Build gene_id -> symbol. Entrez IDs map via NCBI annotation; else id=symbol."""
    idx = counts.index.astype(str)
    if idx.str.fullmatch(r"\d+").mean() > 0.8:  # Entrez gene IDs
        print("[download] gene IDs look like Entrez; downloading annotation for symbols")
        ann = http_get(NCBI_ANNOT_URL)
        if ann is not None:
            annot = pd.read_csv(io.BytesIO(gzip.decompress(ann)), sep="\t", index_col=0)
            annot.index = annot.index.astype(str)
            col = "Symbol" if "Symbol" in annot.columns else annot.columns[0]
            sym = annot[col].reindex(idx)
            return pd.DataFrame({"symbol": sym.values}, index=idx)
    # Ensembl or already-symbol: use the ID itself as the label.
    return pd.DataFrame({"symbol": idx}, index=idx)


def main():
    series_bytes = http_get(SERIES_URL)
    if series_bytes is None:
        sys.exit("Could not download the series matrix at {}".format(SERIES_URL))
    gsm_to_group, gsm_to_title = parse_series(series_bytes)

    # Strategy 1: NCBI uniformly processed raw counts.
    print("\n[download] trying NCBI-generated counts")
    counts = None
    data = http_get(NCBI_COUNTS_URL)
    if data is not None:
        counts = read_counts_bytes(data, "ncbi.tsv.gz")
        print("[download] NCBI counts: {} genes x {} samples".format(*counts.shape))

    # Strategy 2: the series' own supplementary counts file.
    if counts is None:
        print("\n[download] NCBI counts unavailable; scanning supplementary files")
        listing = http_get(SUPPL_URL)
        if listing is None:
            sys.exit("Could not list supplementary files at {}".format(SUPPL_URL))
        fname, all_files = find_counts_file(listing.decode("utf-8", "replace"))
        if fname is None:
            print("[download] no raw-counts file found. Supplementary files are:")
            for n in all_files:
                print("   ", n)
            sys.exit("Point the loader at the right counts file and rerun.")
        print("[download] using supplementary file: {}".format(fname))
        fdata = http_get(SUPPL_URL + fname)
        if fdata is None:
            sys.exit("Failed to download {}".format(SUPPL_URL + fname))
        counts = read_counts_bytes(fdata, fname)
        print("[download] suppl counts: {} genes x {} samples".format(*counts.shape))

    # Map count columns to young/old. Columns may be GSM IDs or sample titles.
    title_to_group = {}
    for gsm, title in gsm_to_title.items():
        if title:
            title_to_group[title] = gsm_to_group.get(gsm)
    col_group = {}
    for col in counts.columns:
        c = str(col)
        if c in gsm_to_group and gsm_to_group[c]:
            col_group[col] = gsm_to_group[c]
        elif c in title_to_group and title_to_group[c]:
            col_group[col] = title_to_group[c]
        else:
            g = classify(c)
            if g:
                col_group[col] = g
    keep_cols = [c for c in counts.columns if c in col_group]
    counts = counts[keep_cols]
    meta = pd.DataFrame(
        {"sample": [str(c) for c in keep_cols],
         "condition": [col_group[c] for c in keep_cols]}
    ).set_index("sample")
    counts.columns = [str(c) for c in counts.columns]

    n_young = int((meta["condition"] == "young").sum())
    n_old = int((meta["condition"] == "old").sum())
    print("\n[download] labeled samples: {} young, {} old".format(n_young, n_old))
    if n_young < 2 or n_old < 2:
        sys.exit("Need at least 2 samples per group; check the parsed labels above.")

    # Coerce to integer counts and apply the low-count filter.
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).round().astype(int)
    keep = counts.sum(axis=1) >= MIN_TOTAL_COUNT
    counts = counts[keep]
    print("[download] kept {} genes after low-count filter".format(counts.shape[0]))

    gene_map = attach_symbols(counts)

    os.makedirs(DATA_DIR, exist_ok=True)
    counts.to_csv(os.path.join(DATA_DIR, "counts.csv"))
    meta.to_csv(os.path.join(DATA_DIR, "meta.csv"))
    gene_map.to_csv(os.path.join(DATA_DIR, "gene_map.csv"))
    print("[download] wrote data/counts.csv, data/meta.csv, data/gene_map.csv")


if __name__ == "__main__":
    main()
