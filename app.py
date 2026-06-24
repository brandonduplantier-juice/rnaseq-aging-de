"""
Interactive explorer for the RNA-seq aging DE results.

Written so a non-specialist can use it: plain-language explainers on every chart,
gene discovery aids so you do not need to know a gene name to start, and real gene
descriptions pulled live from NCBI (via mygene.info) rather than written by hand.

Run locally:
  pip install -r requirements.txt
  python -m streamlit run app.py

Reads results/de_results.csv and results/metrics.json, both committed.
"""

import json
import os

import numpy as np
import pandas as pd
import plotly.express as px
import requests
import streamlit as st

APP_VERSION = "1.2.0"
ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
PADJ_CUTOFF = 0.05


@st.cache_data
def load():
    res = pd.read_csv(os.path.join(RESULTS, "de_results.csv"))
    res = res.dropna(subset=["log2FoldChange", "padj", "baseMean"]).copy()
    res["significant"] = res["padj"] < PADJ_CUTOFF
    res["neglog10padj"] = -np.log10(res["padj"].clip(lower=1e-300))
    res["symbol"] = res["symbol"].astype(str)
    with open(os.path.join(RESULTS, "metrics.json")) as fh:
        meta = json.load(fh)
    return res, meta


@st.cache_data(show_spinner=False)
def gene_info(symbol):
    """Real gene name and summary from NCBI via mygene.info. Never fabricated.

    Returns (full_name, summary, entrez_id). Any field may be None if the lookup
    fails or NCBI has no summary, in which case the UI falls back to link-outs.
    """
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": "symbol:{}".format(symbol), "species": "human",
                    "fields": "name,summary,entrezgene", "size": 1},
            timeout=6,
        )
        r.raise_for_status()
        hits = r.json().get("hits", [])
        if not hits:
            return None, None, None
        h = hits[0]
        return h.get("name"), h.get("summary"), h.get("entrezgene")
    except Exception:
        return None, None, None


def plain_effect(log2fc, padj):
    fold = 2 ** abs(float(log2fc))
    direction = "higher" if log2fc > 0 else "lower"
    sig = "a statistically reliable change" if padj < PADJ_CUTOFF \
        else "not statistically reliable here"
    return ("In older brains, this gene is about {:.1f}x {} than in younger brains "
            "({:+.2f} on a log2 scale). That is {} (adjusted p = {:.2g})."
            ).format(fold, direction, log2fc, sig, padj)


res, meta = load()

st.title("Aging brain gene explorer")
st.caption("Which genes change their activity between young and old human brains. "
           "Based on public RNA-seq data (GSE104704).")

with st.expander("New here? How to read this", expanded=False):
    st.markdown(
        "- **Genes** are instructions; RNA-seq measures how *active* each gene is.\n"
        "- This project compares gene activity in **young vs old human brain** to "
        "find genes that change with age.\n"
        "- **Fold change** is how much louder or quieter a gene gets. A negative "
        "number means quieter (lower) in old age; positive means louder (higher).\n"
        "- **Adjusted p-value** is how sure we are the change is real and not noise. "
        "Smaller is more sure. Below 0.05 is the usual bar for *significant*.\n"
        "- Of all genes tested, only a small set passes that bar. Those are the "
        "interesting ones."
    )

with st.expander("Key terms"):
    st.markdown(
        "- **Gene**: an instruction in your DNA. At any moment a cell only runs some of "
        "them.\n"
        "- **Gene expression**: how active a gene is, how loudly it is being run.\n"
        "- **RNA-seq**: a technology that measures the activity of every gene in a tissue "
        "sample at once.\n"
        "- **Differential expression**: finding genes whose activity differs between two "
        "groups, here young vs old brain, beyond random noise.\n"
        "- **Fold change**: how much louder or quieter a gene gets. Negative means lower "
        "in old age, positive means higher.\n"
        "- **Adjusted p-value (FDR)**: how confident a change is real after testing "
        "thousands of genes at once. Below 0.05 is the usual bar for significant.\n"
        "- **Significant gene**: one that clears that bar, a real age change rather than "
        "chance."
    )

c1, c2, c3, c4 = st.columns(4)
c1.metric("Genes tested", "{:,}".format(meta["n_genes_tested"]),
          help="Genes with enough signal to test for an age difference.")
c2.metric("Significant", meta["n_significant"],
          help="Genes whose change with age is statistically reliable (adjusted p < 0.05).")
c3.metric("Higher in old", meta["n_up_in_old"],
          help="Significant genes that get more active with age.")
c4.metric("Lower in old", meta["n_down_in_old"],
          help="Significant genes that get less active with age. Most age changes here are losses.")

st.subheader("Look up a gene")
st.caption("Do not know a gene name? Pick from the strongest changes below the box, "
           "or just start typing in the box to search all tested genes.")

sig_sorted = res[res["significant"]].sort_values("padj")
top_syms = sig_sorted["symbol"].head(12).tolist()
all_syms = res.sort_values("padj")["symbol"].tolist()
default_idx = all_syms.index(top_syms[0]) if top_syms else 0

choice = st.selectbox("Gene (type to search all tested genes)", all_syms,
                      index=default_idx,
                      help="Start typing to filter. Genes are ordered strongest first.")

if choice:
    row = res[res["symbol"] == choice].iloc[0]
    a, b, c = st.columns(3)
    a.metric("Fold change (log2)", "{:.2f}".format(row["log2FoldChange"]),
             help="Negative = quieter in old age, positive = louder. Each 1.0 doubles the change.")
    b.metric("Adjusted p-value", "{:.2g}".format(row["padj"]),
             help="How sure we are the change is real. Smaller is more sure; under 0.05 is significant.")
    c.metric("Typical activity", "{:.0f}".format(row["baseMean"]),
             help="Average expression level across samples. Higher means a busier gene.")
    st.info(plain_effect(row["log2FoldChange"], row["padj"]))

    name, summary, entrez = gene_info(choice)
    st.markdown("**What this gene does**")
    if name:
        st.write("{} ({})".format(name, choice))
    if summary:
        st.write(summary)
    else:
        st.caption("No plain-language summary was returned for this gene. Use the "
                   "links below for the full annotation.")
    links = []
    if entrez:
        links.append("[NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/{})".format(entrez))
    links.append("[GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene={})".format(choice))
    st.markdown("Read more: " + "  |  ".join(links))
    st.caption("Gene name and summary come from NCBI via mygene.info, not written by us.")

with st.expander("Browse the strongest changes"):
    show = sig_sorted.head(25).copy()
    show["effect"] = show["log2FoldChange"].apply(
        lambda x: "about {:.1f}x {}".format(2 ** abs(x), "higher" if x > 0 else "lower"))
    st.dataframe(
        show[["symbol", "effect", "log2FoldChange", "padj"]],
        use_container_width=True, height=320, hide_index=True,
    )

st.subheader("The big picture (volcano plot)")
st.caption("Each dot is a gene. Left = quieter in old age, right = louder. Higher up "
           "= more statistically reliable. Red dots are the significant ones. Hover any "
           "dot to see the gene.")
fig = px.scatter(
    res, x="log2FoldChange", y="neglog10padj", color="significant",
    color_discrete_map={True: "#cb181d", False: "#bdbdbd"},
    custom_data=["symbol"],
    labels={"log2FoldChange": "quieter in old   <-   fold change   ->   louder in old",
            "neglog10padj": "more reliable (higher)"},
)
fig.update_traces(hovertemplate="<b>%{customdata[0]}</b><br>"
                  "fold change (log2): %{x:.2f}<br>reliability score: %{y:.1f}<extra></extra>")
fig.add_hline(y=-np.log10(PADJ_CUTOFF), line_dash="dash", line_color="black",
              annotation_text="significance line", annotation_position="top left")
fig.update_layout(height=520, legend_title_text="significant")
st.plotly_chart(fig, use_container_width=True)

st.subheader("Full results table")
only_sig = st.checkbox("Show only significant genes", value=True,
                       help="Significant = the change with age is statistically reliable.")
search = st.text_input("Filter by gene name contains", value="").strip().upper()
tbl = res.copy()
if only_sig:
    tbl = tbl[tbl["significant"]]
if search:
    tbl = tbl[tbl["symbol"].str.upper().str.contains(search)]
tbl = tbl.sort_values("padj")[["symbol", "baseMean", "log2FoldChange", "padj", "significant"]]
st.dataframe(tbl, use_container_width=True, height=380, hide_index=True)

st.caption("Population-level research analysis of public data. Not a clinical or "
           "diagnostic tool. App v{}.".format(APP_VERSION))
