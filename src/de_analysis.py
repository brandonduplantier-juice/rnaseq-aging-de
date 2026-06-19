"""
Differential expression: old versus young HSCs, using DESeq2 (via pydeseq2).

DESeq2 models each gene's integer counts with a negative binomial distribution,
estimates how counts differ between groups, shrinks noisy variance estimates by
sharing information across genes, and tests each gene with multiple-testing
correction. It is the standard method for bulk RNA-seq with small sample sizes.

Reads:
  data/counts.csv     genes x samples raw counts
  data/meta.csv       sample, condition (young or old)
  data/gene_map.csv   gene_id, symbol

Writes:
  results/de_results.csv      gene_id, symbol, baseMean, log2FoldChange, lfcSE,
                              stat, pvalue, padj (sorted by padj)
  results/metrics.json        summary counts
  data/normalized_counts.csv  DESeq2-normalized counts (for the plots)

Positive log2FoldChange means higher expression in old than young.
"""

import json
import os

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
RESULTS_DIR = os.path.join(ROOT, "results")
PADJ_CUTOFF = 0.05


def main():
    counts = pd.read_csv(os.path.join(DATA_DIR, "counts.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, "meta.csv"), index_col=0)
    gene_map = pd.read_csv(os.path.join(DATA_DIR, "gene_map.csv"), index_col=0)

    # GEO gene IDs are integers, but pydeseq2 stores feature names as strings.
    # Normalize both sides to str up front so the symbol merge stays aligned.
    counts.index = counts.index.astype(str)
    gene_map.index = gene_map.index.astype(str)

    # DESeq2 wants samples as rows and genes as columns; our matrix is the other
    # way round, so transpose. Align metadata to the same sample order.
    counts_t = counts.T
    meta = meta.loc[counts_t.index]
    print("[de] {} samples, {} genes".format(counts_t.shape[0], counts_t.shape[1]))

    # Fit the DESeq2 model with condition as the only factor.
    dds = DeseqDataSet(counts=counts_t, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()

    # Test old vs young. The contrast order sets the sign: positive log2FC = up
    # in old. summary() runs the Wald test and Benjamini-Hochberg correction.
    stats = DeseqStats(dds, contrast=["condition", "old", "young"], quiet=True)
    stats.summary()
    res = stats.results_df.copy()

    # Attach gene symbols for readability.
    res.insert(0, "symbol", gene_map["symbol"].reindex(res.index))
    res.index.name = "gene_id"
    res = res.sort_values("padj")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    res.to_csv(os.path.join(RESULTS_DIR, "de_results.csv"))

    # Save normalized counts for the plots (samples x genes -> genes x samples).
    normed = pd.DataFrame(
        dds.layers["normed_counts"], index=counts_t.index, columns=counts_t.columns
    ).T
    normed.index.name = "gene_id"
    normed.to_csv(os.path.join(DATA_DIR, "normalized_counts.csv"))

    sig = res[res["padj"] < PADJ_CUTOFF]
    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())
    metrics = {
        "dataset": "GSE104406 (Adelman 2019, human HSC, young vs old)",
        "n_samples": int(counts_t.shape[0]),
        "n_genes_tested": int(res["padj"].notna().sum()),
        "padj_cutoff": PADJ_CUTOFF,
        "n_significant": int(len(sig)),
        "n_up_in_old": n_up,
        "n_down_in_old": n_down,
        "top_genes": sig.head(15)["symbol"].dropna().tolist(),
    }
    with open(os.path.join(RESULTS_DIR, "metrics.json"), "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    print("[de] genes tested: {}".format(metrics["n_genes_tested"]))
    print("[de] significant (padj < {}): {}  ({} up in old, {} down)".format(
        PADJ_CUTOFF, len(sig), n_up, n_down))
    if metrics["top_genes"]:
        print("[de] top genes: {}".format(", ".join(metrics["top_genes"][:10])))
    print("[de] wrote results/de_results.csv and results/metrics.json")


if __name__ == "__main__":
    main()
