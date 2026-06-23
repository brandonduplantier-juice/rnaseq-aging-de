"""
Plots for the young-versus-old DE analysis.

Outputs:
  results/pca.png       sample PCA on normalized counts, colored by group
  results/volcano.png   log2 fold change vs significance, DE genes highlighted
  results/ma_plot.png   mean expression vs fold change, significant genes in red
  results/heatmap.png   z-scored expression of the top DE genes across samples

Reads results/de_results.csv, data/normalized_counts.csv, data/meta.csv.
Run de_analysis.py first.
"""

import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # save files without a display
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
RESULTS_DIR = os.path.join(ROOT, "results")
TOP_N = 40           # genes shown in the heatmap
PADJ_CUTOFF = 0.05

COLORS = {"young": "#2c7fb8", "old": "#d95f0e"}


def main():
    res = pd.read_csv(os.path.join(RESULTS_DIR, "de_results.csv"), index_col=0)
    normed = pd.read_csv(os.path.join(DATA_DIR, "normalized_counts.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, "meta.csv"), index_col=0)
    meta = meta.loc[normed.columns]  # align sample order

    # log2(normalized + 1) is the standard scale for visualizing counts.
    log_expr = np.log2(normed + 1.0)

    # ---- PCA on the most variable genes ----
    top_var = log_expr.var(axis=1).sort_values(ascending=False).head(2000).index
    x = log_expr.loc[top_var].T.values
    x = x - x.mean(axis=0)                       # center genes
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var_explained = (s ** 2 / (s ** 2).sum())[:2] * 100

    fig, ax = plt.subplots(figsize=(6, 5))
    for grp in ("young", "old"):
        idx = [i for i, sm in enumerate(normed.columns) if meta.loc[sm, "condition"] == grp]
        ax.scatter(pcs[idx, 0], pcs[idx, 1], s=60, alpha=0.85,
                   color=COLORS[grp], label=grp, edgecolor="white")
    ax.set_xlabel("PC1 ({:.1f}% var)".format(var_explained[0]))
    ax.set_ylabel("PC2 ({:.1f}% var)".format(var_explained[1]))
    ax.set_title("Sample PCA: young vs old human cortex")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "pca.png"), dpi=150)
    plt.close(fig)

    # ---- Volcano ----
    r = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    r["neglog10padj"] = -np.log10(r["padj"].clip(lower=1e-300))
    sig = r["padj"] < PADJ_CUTOFF
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(r.loc[~sig, "log2FoldChange"], r.loc[~sig, "neglog10padj"],
               s=8, alpha=0.4, color="#bdbdbd", edgecolor="none", label="ns")
    ax.scatter(r.loc[sig, "log2FoldChange"], r.loc[sig, "neglog10padj"],
               s=10, alpha=0.7, color="#cb181d", edgecolor="none",
               label="padj < {}".format(PADJ_CUTOFF))
    ax.axhline(-np.log10(PADJ_CUTOFF), color="black", lw=0.8, ls="--")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("log2 fold change (old vs young)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("Volcano: differential expression with age")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "volcano.png"), dpi=150)
    plt.close(fig)

    # ---- Heatmap of top DE genes ----
    top = res.dropna(subset=["padj"]).sort_values("padj").head(TOP_N).index
    sub = log_expr.loc[top]
    # z-score each gene across samples so colors reflect relative expression.
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    order = [sm for grp in ("young", "old")
             for sm in normed.columns if meta.loc[sm, "condition"] == grp]
    z = z[order]
    labels = res.loc[top, "symbol"].fillna(res.loc[top].index.to_series().astype(str))

    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([meta.loc[sm, "condition"][0] for sm in order], fontsize=7)
    split = sum(1 for sm in order if meta.loc[sm, "condition"] == "young")
    ax.axvline(split - 0.5, color="black", lw=1.5)
    ax.set_title("Top {} DE genes (z-scored, young left, old right)".format(len(top)))
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="z-score")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "heatmap.png"), dpi=150)
    plt.close(fig)

    # ---- MA plot: mean expression vs fold change ----
    m = res.dropna(subset=["padj", "log2FoldChange", "baseMean"]).copy()
    m = m[m["baseMean"] > 0]
    msig = m["padj"] < PADJ_CUTOFF
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(np.log10(m.loc[~msig, "baseMean"]), m.loc[~msig, "log2FoldChange"],
               s=7, alpha=0.35, color="#bdbdbd", edgecolor="none", label="ns")
    ax.scatter(np.log10(m.loc[msig, "baseMean"]), m.loc[msig, "log2FoldChange"],
               s=12, alpha=0.75, color="#cb181d", edgecolor="none",
               label="padj < {}".format(PADJ_CUTOFF))
    ax.axhline(0, color="black", lw=0.8)
    top_ma = m.loc[msig].reindex(
        m.loc[msig, "log2FoldChange"].abs().sort_values(ascending=False).index).head(6)
    for _, row in top_ma.iterrows():
        ax.annotate(str(row["symbol"]),
                    (np.log10(row["baseMean"]), row["log2FoldChange"]),
                    fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("log10 mean expression (baseMean)")
    ax.set_ylabel("log2 fold change (old vs young)")
    ax.set_title("MA plot: expression level vs fold change with age")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "ma_plot.png"), dpi=150)
    plt.close(fig)

    print("[plots] wrote results/pca.png, results/volcano.png, results/ma_plot.png, results/heatmap.png")


if __name__ == "__main__":
    main()
