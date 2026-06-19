# RNA-seq Aging DE: young vs old hematopoietic stem cells

Differential gene expression between young and old human hematopoietic stem cells
(HSCs), using the standard DESeq2 workflow. Built as a portfolio bioinformatics
project: it takes a published aging dataset from raw counts to an annotated list
of differentially expressed genes and the three figures that always accompany a
DE analysis (PCA, volcano, heatmap).

## Dataset

GSE104406 (Adelman et al. 2019, Cancer Discovery). FACS-sorted bone-marrow HSCs,
bulk RNA-seq, 10 young donors (ages 18 to 30) versus 10 aged donors (ages 65 to
75). A balanced two-group design.

Counts come from NCBI's uniformly processed RNA-seq count matrices, so no read
alignment is required. Sample group labels are parsed from the GEO series matrix
at runtime, not hardcoded.

## Method

DESeq2 (via the pydeseq2 library) models each gene's integer counts with a
negative binomial distribution, shares variance information across genes to
stabilize estimates on a small sample, and tests each gene with Benjamini-Hochberg
multiple-testing correction. Positive log2 fold change means higher expression in
old than young. Genes are pre-filtered to those with at least 10 total reads.

## Run it

```
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe run_all.py
```

The first run downloads the counts, annotation, and series metadata, then runs the
analysis and writes the figures. Reruns skip the download.

## Outputs

- results/de_results.csv: every gene with log2 fold change, p value, adjusted p,
  and gene symbol, sorted by significance
- results/metrics.json: summary counts (genes tested, significant, up, down)
- results/pca.png: sample PCA, colored by age group
- results/volcano.png: fold change vs significance
- results/heatmap.png: z-scored expression of the top differentially expressed genes

## Limitations

NCBI-pipeline counts may differ slightly from the authors' own processing, so
exact gene-level numbers can differ from the publication. This is a two-group
analysis with no covariate adjustment (for example sex or batch), which a more
complete study would model. n is 10 per group, typical for this kind of study but
small, so modest effects may not reach significance.
