# Multi-Omics Integration in Colorectal Cancer

Data Science in Life Sciences, Freie Universitat Berlin, 2026
Authors: Ali Raza, Pritom, Niloufer
Supervisor: Kathrina Jahn

## What this is

386 TCGA-COAD/READ patients, three omics layers (RNA-seq, methylation,
mutations), integrated with MOFA2 to find latent factors and check which
ones actually track survival. The signal that comes out of that gets
validated three ways: a Cox model, an independent scRNA-seq cohort, and a
couple of plain machine learning baselines.

Question we're asking?
What molecular programs are associated with colorectal cancer survival, and which cell types drive these programs?
Sub-Questions:
             1-  Which biological processes are associated with poor and favorable prognosis in CRC?
             2-  Which cell types in the tumor  contribute to these prognostic processes?

## Data

| Layer | Source | Patients | Features |
|---|---|---|---|
| RNA-seq | TCGA STAR-Counts, VST-normalized | 647 | 32,313 genes |
| Methylation | Illumina 450K, converted to M-values | 410 | 259,858 CpGs |
| Mutations | Masked somatic MAF, binary matrix | 602 | 100 genes |
| Clinical | TCGA GDC | 633 | survival time, vital status, stage, age |

386 patients have complete data across all three layers after QC and
alignment. None of the raw data is stored in this repo (too large for
GitHub) - it's pulled from these public sources instead:

- TCGA-COAD / TCGA-READ (RNA-seq, methylation, mutations, clinical):
  [GDC Data Portal](https://portal.gdc.cancer.gov/projects/TCGA-COAD),
  [TCGA-READ](https://portal.gdc.cancer.gov/projects/TCGA-READ) - downloaded
  programmatically via `TCGAbiolinks` in `01_QC_EDA.Rmd` (see Reproducing
  below).
- GSE132465 (independent scRNA-seq validation cohort, Lee et al. 2020):
  [NCBI GEO accession GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) -
  download the cell annotation and expression matrix files from there and
  place them in `data/GSE132465/` before running `03_external_scrna_validate.Rmd`.

## Pipeline

1. `01_QC_EDA.Rmd` - download, per-layer QC/normalization, batch
   correction, patient alignment.
2. `02_mofa_integration.Rmd` - MOFA2 (15 factors), k-means subtyping,
   log-rank survival screening per factor.
3. `03_external_scrna_validate.Rmd` - do the top Factor 6/12 genes show
   cell-type-specific expression in an independent scRNA-seq CRC cohort
   (GSE132465)?
4. `04_survival_analysis.Rmd` - Cox model on those genes, Kaplan-Meier
   curves.
5. `05_ml_classification.py` - Random Forest / XGBoost predicting 3-year
   survival directly, with class-weighting, SMOTE, and undersampling
   compared side by side.
6. `06_random_survival_forest.py` - same genes, a Random Survival Forest
   instead of Cox.
7. `07_annotation.Rmd` - gene symbols/names for the top Cox genes,
   g:Profiler pathway enrichment.

## Results

- Factor 6 & 12 is the only MOFA2 factor tied to survival.
- Cox model on the Factor 6/12 genes: C-index 0.732, KM p < 0.0001 for
  high vs. low risk.
- 47 of 49 genes detected in the scRNA-seq cohort, with clear cell-type
  specific expression.
- Pathway enrichment points to B-cell mediated immunity.
- Random Forest / XGBoost classifiers (3-year survival, direct
  classification): best around 0.60 accuracy - a modest result, not a
  strong signal on its own.
- Random Survival Forest: C-index 0.573 - similarly modest.

Worth being honest about that last two: the Cox model does noticeably
better than either ML baseline here, most likely because 168-370 patients
and 50 genes just isn't much for a forest model to learn from.

## Repository structure

```
Analysis/     Rmd files and Python scripts listed above
data/         final matrices/tables each step depends on (rds, csv)
Outputs/      one subfolder of figures per analysis file
```

Large raw/intermediate files (`.rds`, GDC/GEO downloads) are excluded via
`.gitignore` and regenerated locally when the pipeline runs.

## Requirements

R 4.x:
```r
BiocManager::install(c("MOFA2", "DESeq2", "TCGAbiolinks",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19", "org.Hs.eg.db",
  "maftools", "sva"))
install.packages(c("survival", "survminer", "cluster", "umap",
  "gprofiler2", "ggplot2", "dplyr", "tidyr", "reshape2", "data.table",
  "Seurat"))
```

Python 3.x:
```
pandas, scikit-learn, xgboost, imbalanced-learn, scikit-survival,
pyreadr, matplotlib
```

## Reproducing

1. Clone the repo.
2. Run `01_QC_EDA.Rmd` first (its download chunks are `eval=FALSE` by
   default - flip to `TRUE` for a first run to pull data from GDC).
3. Knit the rest of the files in order, 02 through 07. MOFA2 training
   (`02_mofa_integration.Rmd`) takes 30-60 minutes and is also gated
   behind `eval=FALSE`; it reloads a cached model after that.

`03_external_scrna_validate.Rmd` expects the GSE132465 files in
`data/GSE132465/`. The two Python scripts (05, 06) need
`data/top_extracted50_genes.csv` and `data/clin_final.rds`, both produced
earlier in the pipeline.

## References

1. Guinney J, et al. The consensus molecular subtypes of colorectal cancer.
   Nature Medicine, 21, 1350-1356 (2015).
2. Argelaguet R, et al. Multi-Omics Factor Analysis. Molecular Systems
   Biology, 14, e8124 (2018).
3. Lee HO, et al. Lineage-dependent gene expression programs influence the
   immune landscape of colorectal cancer. Nature Genetics, 52, 594-603
   (2020).
4. Johnson WE, et al. Adjusting batch effects in microarray expression data
   using empirical Bayes methods. Biostatistics, 8, 118-127 (2007).
