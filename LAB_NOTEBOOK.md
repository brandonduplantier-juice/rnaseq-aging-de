# RNA-seq Aging DE: Lab Notebook and Learning Log

A living document. It explains what this project is, how each piece works, why we
made each choice, and where we are. We bump the version and add a changelog row
every time we change it, then commit. Git stores the real diffs, this header keeps
it readable.

Version: v0.1
Last updated: 2026-06-19
Owner: Brandon

## Version history

| Version | Date       | Change                                            |
|---------|------------|---------------------------------------------------|
| v0.1    | 2026-06-19 | First notebook. Project scaffolded and tested on synthetic data, not yet run on the real dataset. |

How to update this file: make your edits, bump the version number above, add one
changelog row, then commit.

## 1. North Star

The long-term goal is a biomedical longevity company. The aging clock measured
how old a tissue looks. This project asks a different question: which genes
actually change as a tissue ages. Differential expression is the most common
day-to-day task in a bioinformatics analyst role, so this piece covers the
skill the clock did not, while staying on the aging theme by comparing young and
old blood stem cells.

## 2. What this project does, in plain language

It takes a published dataset of blood stem cells from 10 young and 10 old people,
counts how active each gene is in each person, and finds the genes that are
reliably more or less active in the old group. It then draws the three standard
pictures: a PCA that shows whether young and old samples separate, a volcano plot
that shows which genes changed and how strongly, and a heatmap of the top genes.

## 3. Concepts you should know

- Counts: for each gene in each sample, how many sequencing reads landed on it.
  Higher counts mean the gene was more active. The input is a genes-by-samples
  table of these counts.
- Normalization: samples are sequenced to different depths, so raw counts are not
  directly comparable. DESeq2 computes size factors to put samples on a common
  scale before comparing.
- log2 fold change: the difference between groups, on a log2 scale. +1 means
  twice as high in old, -1 means half. We set it so positive is up in old.
- p value and adjusted p (padj): the p value is the chance of seeing a difference
  this large if the gene truly did not change. Because we test thousands of genes,
  we adjust with Benjamini-Hochberg to control false discoveries. We call genes
  significant at padj < 0.05.
- DESeq2: the standard method for this. It models counts with a negative binomial
  distribution and borrows information across genes to get stable estimates when
  the sample size is small. We use pydeseq2, the Python implementation.
- PCA, volcano, heatmap: the three figures that accompany almost every DE study.

## 4. Dataset

GSE104406 (Adelman et al. 2019, Cancer Discovery). FACS-sorted bone-marrow
hematopoietic stem cells, bulk RNA-seq, 10 young (18 to 30) vs 10 aged (65 to 75).
Counts are NCBI's uniformly processed gene-level counts, so no alignment is
needed. Group labels are parsed from the GEO series matrix at runtime.

## 5. Files

- src/download_data.py: downloads the NCBI counts matrix, the gene annotation, and
  the series metadata; parses young vs old labels; writes counts.csv, meta.csv,
  gene_map.csv.
- src/de_analysis.py: runs DESeq2 (old vs young), writes de_results.csv,
  metrics.json, and normalized counts for the plots.
- src/plots.py: writes pca.png, volcano.png, heatmap.png.
- run_all.py: runs all three in order; skips the download if counts.csv exists.

## 6. How to run

```
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe run_all.py
```

This project needs its own virtual environment. pydeseq2 currently requires pandas
below 3.0, while the aging clock used pandas 3.0, so keeping the environments
separate avoids a conflict.

## 7. Limitations

NCBI-pipeline counts can differ slightly from the authors' own processing, so
exact numbers may not match the paper. It is a two-group comparison with no
covariate adjustment (sex, batch), and n is 10 per group, which is typical but
small.

## 8. Current status

Scaffold built. The DE and plotting steps were tested end to end on synthetic
counts and behaved correctly (recovered injected differentially expressed genes,
produced all three figures). Not yet run on the real dataset. Next action: run the
pipeline and record the first numbers here.

Watch items on first run: the download is from NCBI and may be slow; the group
labels are auto-parsed, so check the printed "labeled samples: N young, N old"
line reads 10 and 10.

## 9. Open questions

- How many genes are significant at padj < 0.05, and how many up vs down in old.
- Do young and old samples separate cleanly in the PCA.
- Do known HSC-aging genes (for example EGR1, reported up in old in this dataset)
  appear among the top hits.

## 10. Next steps

- Run on the real data and record metrics here (v0.2).
- Add a gene-set enrichment step on the significant genes (which pathways change
  with age) as a strong follow-up.
- Push to GitHub and add a portfolio card.
