# RNA-seq Aging Differential Expression: Young vs Old Human Brain

Finds which genes change their activity level between young and old human cortex,
using the standard DESeq2 workflow on public data, fully reproducible.

Version 1.2.0

## The result

**119 genes** are differentially expressed at 5 percent FDR between young and old
human cortex (GSE104704, 8 young and 10 aged donors). The shift is strongly
one-sided: **113 genes go down with age and only 6 go up.** Several of the strongest
are myelin and oligodendrocyte genes (ZNF488, GJC2, CNDP1), which lines up with the
well-documented decline of myelin in the aging brain. A pathway-enrichment pass finds
9 biological-process terms significant at 5 percent FDR.

![MA plot](results/ma_plot.png)

## In plain English

No biology background needed. Here is the idea, with each term defined as it comes up.

**What a gene "doing" something means.** Every cell carries the same DNA, but at any
moment only some genes are switched on. *Gene expression* is how active a gene is,
how much its instructions are being read and turned into working molecules. Two cells
with identical DNA can behave very differently because different genes are turned up
or down.

**What RNA-seq measures.** When a gene is active it produces RNA copies. *RNA-seq* is
a technology that counts those copies for every gene at once, giving a number per
gene that reflects how active it was. More copies means more activity.

**What this project asks.** Take brain tissue from young donors and old donors, and
find the genes whose activity reliably differs between the two groups. A gene that
differs is called *differentially expressed*. Those genes are candidates for what
changes in the brain as it ages.

**Why you cannot just compare averages.** With only 8 and 10 donors and thousands of
genes, random noise will make plenty of genes look different by luck. The analysis
uses *DESeq2*, a statistical method built for exactly this, and applies a *multiple-
testing correction* (FDR, the false discovery rate) so that "significant" means
unlikely to be a fluke, not just different on paper.

**What the result says.** 119 genes pass that bar, and almost all of them go *down*
with age. The strongest ones relate to myelin, the insulating sheath around nerve
fibers, whose breakdown in aging is already well established. So the analysis
recovered a known biological story from raw data, which is the point.

What this is **not**: it is a computational analysis of public data, not wet-lab
work, and a gene changing with age is an association, not proof it causes aging.

## How it works (technical)

DESeq2 (via pydeseq2) models each gene's integer read counts with a negative
binomial distribution, shares variance information across genes to stabilize
estimates on small samples, and applies Benjamini-Hochberg FDR correction. Genes are
pre-filtered to at least 10 total reads. A positive log2 fold change means higher in
old. Pathway enrichment runs against Enrichr gene-set libraries.

## What is in this repo

    run_all.py           one command: download, analyze, write table and figures
    app.py               interactive gene explorer (Streamlit)
    src/                 loader, DESeq2 wrapper, enrichment, plotting
    results/             de_results.csv, metrics.json, ma_plot.png, volcano, enrichment
    data/                empty by design; first run downloads the counts

## Explore it

    pip install -r requirements.txt
    streamlit run app.py

Search any gene to see its age-related change in plain language, hover an interactive
volcano plot, and filter the full results table.

## Reproduce

    python -m venv .venv
    .venv\Scripts\python.exe -m pip install -r requirements.txt
    .venv\Scripts\python.exe run_all.py

First run downloads the counts, then writes the table, metrics, and figures. Reruns
skip the download. Pathway enrichment needs internet (Enrichr):

    .venv\Scripts\python.exe run_all.py --enrichment

## Limitations

Small cohort (18 donors), one brain region, bulk tissue (it measures the average
across many cell types, not individual cells). Differential expression shows
association with age, not causation. Gene-family interpretations are framed as
consistent with known brain aging, not as discoveries.

## Glossary

- **Gene expression**: how active a gene is, how much it is being read and used.
- **RNA-seq**: a method that counts the RNA copies a gene produces, as a readout of
  activity.
- **Differentially expressed**: reliably more or less active in one group than
  another.
- **DESeq2**: the standard statistical tool for finding those genes from count data.
- **Counts**: the raw number of RNA reads measured per gene.
- **FDR (false discovery rate)**: a correction so "significant" results are unlikely
  to be flukes; 5 percent FDR means about 5 percent of called genes may be false.
- **Log2 fold change**: the size and direction of the change; positive means higher
  in old.
- **Volcano / MA plot**: standard charts that show effect size against significance.

## Citation and disclaimer

Data: GSE104704 (Nativio et al., human brain). Method: Love et al., DESeq2, Genome
Biol 2014. Research and portfolio project.
