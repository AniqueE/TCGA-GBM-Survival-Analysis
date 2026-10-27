# Transcriptomic analysis of age and survival in TCGA-GBM

Code and processed data for an honours research project in Bioinformatics at the
University of Cape Town.

The study stratifies glioblastoma patients from The Cancer Genome Atlas (TCGA) by
age at diagnosis and survival time, identifies differentially expressed genes
between the resulting groups, characterises the affected biological processes, and
derives a candidate prognostic gene signature.

**Author:** Anique Engelbrecht
**Supervisor:** Dr Hocine Bendou
**Co-supervisor:** Dr Michelle Livesey
**Institution:** University of Cape Town

---

## Data

Bulk RNA-sequencing expression values and clinical phenotype data for TCGA-GBM were
obtained from the [UCSC Xena platform](https://xenabrowser.net/) (Goldman et al., 2020).

Patients were retained if they had no prior cancer diagnosis or treatment, were
deceased at the time of data collection, and had a primary tumour with tissue origin
recorded as "Brain, NOS" and sample type classified as tumour. One patient
(TCGA-06-0211) contributed two vials from the same tumour; only the first (01A) was
kept so that each patient appears once.

This yielded **122 patients** and **19,414 protein-coding genes**.

Patients were assigned to groups using an age cutoff of 55 years and a survival
cutoff of 365 days:

|           | Long survival (LS) | Short survival (SS) |
|-----------|--------------------|---------------------|
| **Old**   | 40                 | 45                  |
| **Young** | 22                 | 15                  |

The filtered expression matrix and phenotype table used in this study are in `data/`.
Both are TCGA open-access data; use is subject to the
[TCGA data use policies](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing).

---

## Repository structure

```
GBM_analysis.R   the full analysis pipeline, run top to bottom
data/            filtered expression matrix and phenotype table (.csv.gz)
results/         output tables (DEG lists, enrichment results, survival metrics)
figures/         figures as they appear in the write-up
```

---

## Running the analysis

The whole pipeline is in a single script, `GBM_analysis.R`, which is run from top to
bottom in a clean R session. Each section depends on objects created by the sections
before it, so the script cannot be run out of order. Sections are marked with banner
comments and appear in this order:

| # | Section |
|---|---------|
| 1 | Packages and setup |
| 2 | Data acquisition, filtering and group assignment |
| 3 | Gene symbol and biotype annotation (biomaRt) |
| 4 | Count back-transformation and DESeq2 dataset construction |
| 5 | Differential gene expression and the four pairwise comparisons |
| 6 | Volcano plots |
| 7 | Venn diagrams |
| 8 | GO and KEGG over-representation analysis |
| 9 | Hallmark gene set enrichment analysis |
| 10 | Variance-stabilising transformation and candidate gene set |
| 11 | Univariate Cox proportional hazards regression |
| 12 | LASSO-Cox regression |
| 13 | Risk score model and evaluation plots |s
| 14 | Signature gene Kaplan-Meier curves |
| 15 | Signature gene boxplots |
| 16 | Multivariate Cox model and proportional hazards assumption |
| 17 | Internal bootstrap validation |

### Analysis summary

**Differential expression.** DESeq2 with an interaction design, Wald test, no log2
fold change shrinkage, independent filtering disabled, and Benjamini-Hochberg
adjustment. Genes were considered significant at padj ≤ 0.05 and |log2FC| ≥ 1,
giving 156 (Old: SS vs LS), 23 (Young: SS vs LS), 15 (SS: Young vs Old) and
73 (LS: Young vs Old) differentially expressed genes.

**Functional enrichment.** Over-representation analysis against GO (all three
ontologies) and KEGG using clusterProfiler, with the 19,414 protein-coding genes as
the background universe. Gene set enrichment analysis against the 50 MSigDB Hallmark
gene sets, ranking genes by the DESeq2 Wald statistic.

**Survival analysis.** The union of the four differentially expressed gene lists
(237 genes) was used as input to univariate Cox regression and LASSO-Cox regression.
LASSO-Cox at lambda.min retained six genes (PTPRN, CLDN16, IGFBP6, KHDRBS2, OR10Q1
and SARNP). A risk score was computed as the weighted sum of variance-stabilising
transformed expression, and evaluated with Kaplan-Meier analysis, time-dependent ROC
curves, a multivariate Cox model adjusted for age, and internal bootstrap validation.

---

## Reproducibility

A seed of **123** is set before LASSO-Cox model fitting and before the bootstrap
procedure. The seed option is enabled for gene set enrichment analysis, which uses
permutation-based p-value estimation and would otherwise vary between runs.

### Software

Analyses were run in R version 4.5.2.

| Package | Version |
|---------|---------|
| DESeq2 | 1.50.2 |
| biomaRt | 2.66.1 |
| EnhancedVolcano | 1.28.2 |
| ggVennDiagram | 1.5.7 |
| clusterProfiler | 4.18.4 |
| enrichplot | 1.30.5 |
| org.Hs.eg.db | 3.22.0 |
| msigdbr | 26.1.0 |
| survival | 3.8-9 |
| survminer | 0.5.2 |
| glmnet | 5.0 |
| timeROC | 0.4.1 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |

Full session details can be regenerated with `sessionInfo()`.

---

## Citation

If you use this code, please cite the associated thesis:

> Anique Engelbrecht (2026). *Glioblastoma Survival Analysis Considering Age as a Cofactor and Using Bulk RNA-Sequencing Data*. Honours thesis, University of Cape Town.

The underlying data should be cited as:

> Goldman, M. J., Craft, B., Hastie, M., Repečka, K., McDade, F., Kamath, A., ... &
> Haussler, D. (2020). Visualizing and interpreting cancer genomics data via the Xena
> platform. *Nature Biotechnology*, *38*(6), 675–678.

---

## Licence

Code in this repository is released under the MIT Licence (see `LICENSE`).
The TCGA data in `data/` remain subject to the TCGA data use policies linked above.
