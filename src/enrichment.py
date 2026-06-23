"""
Pathway enrichment for the significant DE genes.

Runs over-representation analysis (Enrichr via gseapy) on the genes that pass
padj < 0.05 in results/de_results.csv, then writes:
  results/enrichment.csv   ranked terms with adjusted p values
  results/enrichment.png   bar chart of the top terms by -log10 adjusted p

Needs network access to Enrichr, so run on a machine with internet:
  python src/enrichment.py

Honest-null behavior: if no term passes the FDR cutoff, the script says so and
still plots the top suggestive terms, clearly labeled as not significant.
Reporting a null is a real result, not a failure.
"""

import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
PADJ_CUTOFF = 0.05
FDR_CUTOFF = 0.05
GENE_SETS = ["GO_Biological_Process_2021"]
TOP_N = 15


def main():
    res = pd.read_csv(os.path.join(RESULTS, "de_results.csv"))
    sig = res.dropna(subset=["padj", "symbol"])
    sig = sig[sig["padj"] < PADJ_CUTOFF]
    genes = sorted(set(sig["symbol"].astype(str)))
    if len(genes) < 5:
        print("[enrichment] too few significant genes to test: {}".format(len(genes)))
        return
    print("[enrichment] testing {} significant genes against {}".format(
        len(genes), ", ".join(GENE_SETS)))

    try:
        import gseapy as gp
    except ImportError:
        print("[enrichment] gseapy not installed. Run: pip install gseapy")
        return

    try:
        enr = gp.enrichr(gene_list=genes, gene_sets=GENE_SETS, outdir=None)
    except Exception as exc:
        print("[enrichment] Enrichr call failed (needs internet): {}".format(exc))
        return

    df = enr.results.sort_values("Adjusted P-value").reset_index(drop=True)
    df.to_csv(os.path.join(RESULTS, "enrichment.csv"), index=False)

    n_sig = int((df["Adjusted P-value"] < FDR_CUTOFF).sum())
    if n_sig == 0:
        print("[enrichment] honest null: 0 terms pass FDR < {}. "
              "Plotting top suggestive terms, labeled not significant.".format(FDR_CUTOFF))
    else:
        print("[enrichment] {} terms pass FDR < {}".format(n_sig, FDR_CUTOFF))

    top = df.head(TOP_N).iloc[::-1]
    vals = -np.log10(top["Adjusted P-value"].clip(lower=1e-300))
    colors = ["#cb181d" if p < FDR_CUTOFF else "#bdbdbd"
              for p in top["Adjusted P-value"]]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(range(len(top)), vals, color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([str(t)[:55] for t in top["Term"]], fontsize=7)
    ax.axvline(-np.log10(FDR_CUTOFF), color="black", lw=0.8, ls="--")
    ax.set_xlabel("-log10 adjusted p")
    title = "Pathway enrichment of age-DE genes"
    if n_sig == 0:
        title = title + " (no term significant at FDR 0.05)"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "enrichment.png"), dpi=150)
    plt.close(fig)
    print("[enrichment] wrote results/enrichment.csv and results/enrichment.png")


if __name__ == "__main__":
    main()
