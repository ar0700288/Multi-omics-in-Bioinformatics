# 🧬 Multi-Omics Integration in Colorectal Cancer
### Survival Prediction via MOFA2 · Cox Regression · scRNA-seq Validation

<p align="center">
  <img src="https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white"/>
  <img src="https://img.shields.io/badge/Framework-MOFA2-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Data-TCGA--COAD%2FREAD-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Status-Complete-brightgreen?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-MIT-blue?style=for-the-badge"/>
</p>

---

> **Course:** Data Science in Life Sciences · Freie Universität Berlin · 22026  
> **Authors:** Ali Raza · Pritom · niloufer 
> **Supervisor:** [Kathrina Jahn, FU Berlin]

---

## 📌 Overview

This project implements a complete **multi-omics integration pipeline** for colorectal cancer (CRC) using data from **386 TCGA-COAD/READ patients**. By jointly integrating RNA-seq, DNA methylation, and somatic mutation data through **MOFA2**, we identify survival-associated latent factors and validate them through Cox regression and independent single-cell RNA-seq analysis.

### 🎯 Research Question
> *Can multi-omics integration reveal survival-associated molecular signatures in CRC that single-layer analysis cannot detect?*

---

## 🔬 Key Results

| Metric | Value |
|--------|-------|
| Patients (aligned) | 386 |
| Survival-significant factor | **Factor 12 (p = 0.020)** |
| Cox C-index | **0.732** |
| KM High vs Low risk | **p < 0.0001** |
| scRNA-seq genes validated | **47 / 49** |
| Top pathway | **B cell mediated immunity** |

---

## 🗂️ Repository Structure

```
crc-multiomics-integration/
│
├── 📄 README.md                          ← You are here
├── 📄 .gitignore                         ← Excludes large RDS files
│
├── 📓 analysis/
│   └── crc_multiomics.Rmd               ← Main analysis pipeline
│
├── 📊 data/
│   ├── top_genes_Factor6.csv            ← Top 50 survival genes (Factor 12 + 6)
│   ├── ml_f6_features.csv               ← ML feature matrix (386 patients × 50 genes)
│   ├── top15_annotated_genes.csv        ← Annotated Cox survival genes
│   └── gprofiler_results.csv            ← Pathway enrichment results
│
└── 🖼️ outputs/
    ├── factor6_heatmap.png              ← scRNA-seq gene expression heatmap
    ├── gprofiler_barplot.png            ← Pathway enrichment bar plot
    ├── rsf_importance.png               ← Random Survival Forest importance
    ├── oncoprint.png                    ← Mutation landscape OncoPrint
    ├── pathway_table.png                ← Enriched pathways table
    └── annotation_table.png             ← Annotated survival genes table
```

---

## 🔄 Analysis Pipeline

```
TCGA-COAD + TCGA-READ
         │
         ▼
┌─────────────────────────────────────┐
│  Phase 1 — Data Acquisition         │
│  TCGAbiolinks: RNA + Meth + Mut     │
└─────────────────┬───────────────────┘
                  │
                  ▼
┌─────────────────────────────────────┐
│  Phase 2 — Preprocessing & QC       │
│  VST · M-values · Binary mutations  │
│  ComBat batch correction            │
│  Patient alignment (n=386)          │
└─────────────────┬───────────────────┘
                  │
                  ▼
┌─────────────────────────────────────┐
│  Phase 3 — EDA                      │
│  PCA · UMAP · OncoPrint · TMB       │
└─────────────────┬───────────────────┘
                  │
                  ▼
┌─────────────────────────────────────┐
│  Phase 4 — MOFA2 Integration        │
│  15 latent factors                  │
│  RNA (Gaussian) · Meth (Gaussian)   │
│  Mutations (Bernoulli)              │
└────────────┬────────────────────────┘
             │
      ┌──────┴──────┐
      ▼             ▼
 Branch A        Branch B
 k-means         Factor survival
 subtyping       screening
 (k=2)           Factor12 p=0.020 ★
      │             │
      └──────┬───────┘
             ▼
┌─────────────────────────────────────┐
│  Phase 5 — Survival Analysis        │
│  Top 25 genes × 2 factors = 50      │
│  Cox Regression  C-index = 0.732    │
│  RSF             C-index = 0.587    │
│  KM p < 0.0001                      │
└─────────────────┬───────────────────┘
                  │
                  ▼
┌─────────────────────────────────────┐
│  Phase 6 — scRNA-seq Validation     │
│  GSE132465 · 2,500 tumor cells      │
│  47/49 genes validated              │
│  B cell + Epithelial enrichment     │
└─────────────────┬───────────────────┘
                  │
                  ▼
┌─────────────────────────────────────┐
│  Phase 7 — Biological Annotation    │
│  org.Hs.eg.db · g:Profiler          │
│  B cell mediated immunity (p=0.001) │
└─────────────────────────────────────┘
```

---

## 🧪 Methods Summary

### Data
| Layer | Source | Samples | Features |
|-------|--------|---------|----------|
| RNA-seq | TCGA STAR-Counts → VST | 647 | 32,313 genes |
| Methylation | Illumina 450K → M-values | 410 | 259,858 CpGs |
| Mutations | Masked MAF → binary | 602 | 100 genes |
| Clinical | TCGA GDC | 633 | Survival, stage, age |

### Key Methodological Choices
- **Batch correction:** ComBat (sva) — reduced RNA R² from 0.235 → 0.164
- **MOFA2:** Gaussian likelihood for RNA/methylation, Bernoulli for mutations
- **Feature selection:** Top 25 genes per factor by absolute MOFA2 weight
- **Survival model:** Cox PH regression + Random Survival Forest (comparison)
- **External validation:** GSE132465 Korean CRC cohort (Lee et al. 2020)

---

## 📦 Dependencies

```r
# Bioconductor
BiocManager::install(c(
  "MOFA2", "DESeq2", "TCGAbiolinks",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "org.Hs.eg.db", "maftools"
))

# CRAN
install.packages(c(
  "survival", "survminer", "randomForestSRC",
  "gbm", "sva", "cluster", "umap",
  "gprofiler2", "ggplot2", "dplyr",
  "tidyr", "kableExtra", "gridExtra"
))

# Seurat
install.packages("Seurat")
```

---

## 🚀 How to Reproduce

1. **Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/crc-multiomics-integration.git
cd crc-multiomics-integration
```

2. **Download raw data** — Run `eval=FALSE` chunks in `crc_multiomics.Rmd` once to download from TCGA GDC

3. **Run the pipeline** — Knit `analysis/crc_multiomics.Rmd` in RStudio

> ⚠️ Note: MOFA2 training takes 30–60 minutes. The trained model is not included due to file size limits. Set `eval=TRUE` on the `mofa-run` chunk for first run.

---

## 📈 Key Figures

### Factor Survival Association
Factor 12 (p = 0.020) is the only MOFA2 factor significantly associated with survival after ComBat batch correction.

### Cox Survival Stratification  
KM p < 0.0001 — High risk patients drop to 30% survival while Low risk patients maintain >70% at 4000 days.

### scRNA-seq Validation
47/49 Factor genes confirmed in independent Korean CRC cohort. B cell immunoglobulin genes and epithelial markers show cell-type specific expression.

### Pathway Enrichment
B cell mediated immunity and immunoglobulin-mediated immune response are the top enriched pathways (FDR p < 0.01).

---

## 📚 References

1. Guinney J, et al. (2015). The consensus molecular subtypes of colorectal cancer. *Nature Medicine*, 21, 1350–1356.
2. Argelaguet R, et al. (2018). Multi-Omics Factor Analysis. *Molecular Systems Biology*, 14, e8124.
3. Lee HO, et al. (2020). Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. *Nature Genetics*, 52, 594–603.
4. Johnson WE, et al. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8, 118–127.

---

## 📄 License

This project is licensed under the MIT License.

---

<p align="center">
  <i>Freie Universität Berlin · M.Sc. Bioinformatics · Data Science in Life Sciences · 2026</i>
</p>
