# Alzheimer's disease genetic risk in single brain cells

## Summary

**Main finding.** Inherited Alzheimer's disease (AD) risk is concentrated in **microglia**, the
brain's immune cells. The specific microglial and astrocyte states that Green et al. (2024) proposed
as causal drivers of AD (Mic.12, Mic.13 and Ast.10) do not carry clear extra genetic risk beyond
microglia as a whole.

**The paper this project builds on.** Green et al. (2024, *Nature*) studied 1.65 million brain
nuclei from 437 people in the ROSMAP cohort. Using gene expression alone, they named three glial
states as causal drivers of AD: two microglial states (Mic.12, Mic.13) and one astrocyte state
(Ast.10). They did not test whether these states carry inherited genetic risk.

**The question.** Do the cell states Green et al. called causal also carry inherited AD genetic
risk? If yes, genetics would independently support their claim. If no, the states are more likely
a reaction to disease than its inherited starting point.

**Why this dataset.** Green et al.'s single-cell data (ROSMAP) requires an access request through
Synapse / the AD Knowledge Portal, which has not been granted yet. This project therefore uses a
freely available dataset of the same kind as a pilot: GSE138852 (Grubman et al. 2019), 13,214
nuclei from the entorhinal cortex of 6 AD and 6 control donors. The code is written so the same
pipeline can be run on ROSMAP once access is granted.

**What we found.**

1. AD genetic risk genes are strongly active in microglia (group test p = 0.001, z = 6.8) and in no
   other major cell type. The result is the same with the authors' labels, our own clustering and
   the Allen Institute reference labels.
2. The risk is equally high in normal (homeostatic) and disease-associated microglia, and in AD and
   control donors. Risk marks microglia as a cell type, not a disease state.
3. Cells scoring high for Mic.13 first appeared to carry more risk (p = 2e-14 across all nuclei).
   This was mostly because they include more microglia, and because the gene *APOE* is counted on
   both sides of the test. Inside microglia, with *APOE* removed, the correlation disappears
   (rho = 0.03, p = 0.37); only a small top-quartile difference remains (p = 0.008).
4. Mic.12 (without *APOE*) and Ast.10 show no link to genetic risk.

**Conclusion.** Inherited AD risk points to microglia in general. The states Green et al. called
causal look more like reactions to disease than genetically primed starting points. This is a
12-donor pilot with marker-gene proxies for the states, so it needs to be repeated on ROSMAP or
SEA-AD with real state labels.

## Data used

| Data | Source | What it is |
|---|---|---|
| Single-nucleus RNA-seq | GEO [GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852), Grubman et al. 2019 | 13,214 nuclei, 10,850 genes, entorhinal cortex, 6 AD + 6 control |
| AD GWAS | GWAS Catalog [GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158), Bellenguez et al. 2022 | 21.1 million variants, 111,326 cases and 677,663 controls |
| Reference cell labels | Allen Institute [MapMyCells](https://knowledge.brain-map.org/mapmycells/process), SEA-AD MTG taxonomy | Cell type and cell state for every nucleus |
| LD reference | 1000 Genomes European panel (from MAGMA) | Used to correct for correlated variants |

Download instructions are in [data/README.md](data/README.md) and [tools/README.md](tools/README.md).

## What each notebook does

| Notebook | Input | What happens | Main output |
|---|---|---|---|
| [01_preprocessing](notebooks/01_preprocessing.ipynb) | GEO counts + metadata | QC metrics, normalise, log-transform | `interim/adata_lognorm.h5ad` |
| [02_clustering](notebooks/02_clustering.ipynb) | output of 01 | 3,000 variable genes, PCA, UMAP, Leiden clusters | `interim/adata_clustered.h5ad` |
| [03_gwas_qc](notebooks/03_gwas_qc.ipynb) | GWAS file | Remove ambiguous, rare and invalid variants (7.36 M kept) | `interim/magma_input.txt` |
| [04_magma_gene_analysis](notebooks/04_magma_gene_analysis.ipynb) | output of 03 | MAGMA turns variant p-values into one risk Z score per gene | `results/tables/04_magma_genes.tsv` |
| [05_marker_scoring](notebooks/05_marker_scoring.ipynb) | output of 02 | Mic.12, Mic.13, Ast.10 marker-gene scores for every nucleus | `interim/marker_scores.tsv` |
| [06_scdrs](notebooks/06_scdrs.ipynb) | outputs of 02 and 04 | scDRS: AD risk score per nucleus vs 1,000 random gene sets | `interim/scdrs_full_score.tsv` |
| [07_mapmycells_integration](notebooks/07_mapmycells_integration.ipynb) | outputs of 02, 05, 06 + Allen labels | Check labels, test which cell types carry risk, microglial states | `results/tables/07_*.tsv` |
| [08_enrichment_test](notebooks/08_enrichment_test.ipynb) | output of 07 | Test whether Mic.12, Mic.13, Ast.10 carry extra risk | `results/tables/08_state_enrichment.tsv` |

## Results in more detail

### GWAS and gene-level risk

After quality control, 7,356,090 variants remain; 4,683 are genome-wide significant, led by the
*APOE* region on chromosome 19 and *BIN1* on chromosome 2. MAGMA finds 147 genes past the Bonferroni
threshold (p < 2.7e-6), including *APOC1*, *INPP5D*, *MS4A6A*, *BIN1*, *APOE* and *TREM2*.

### Reference labels agree with the authors

The Allen reference labels agree with the authors' cell types for 98.2 % of nuclei (96.8 to 99.8 %
per cell type). They also give confident labels to the 925 "unID" and 405 "doublet" nuclei the
authors could not assign.

### Which cell types carry AD risk

scDRS group test: the top 5 % of scores in each cell type compared with 1,000 random gene sets.

| Cell type (Allen subclass) | Nuclei | Mean risk score | p | z |
|---|---:|---:|---:|---:|
| Microglia | 668 | +1.34 | 0.001 | 6.80 |
| Lamp5 neurons | 129 | +0.30 | 0.007 | 2.86 |
| Oligodendrocytes | 7,514 | +0.03 | 0.24 | 0.71 |
| Astrocytes | 2,461 | -0.18 | 0.66 | -0.44 |
| OPC | 1,394 | -0.53 | 0.97 | -1.84 |

0.001 is the smallest p-value possible with 1,000 random gene sets. Microglia average +1.34 against
-0.07 for all other nuclei.

### Microglial states

| Microglial state | AD | Control | More common in AD? | Mean risk score |
|---|---:|---:|---|---:|
| Micro-PVM_3-SEAAD (disease-associated) | 85 | 34 | yes, 3.3 times (FDR 6e-8) | 1.73 |
| Micro-PVM_2_3-SEAAD (disease-associated) | 118 | 73 | yes, 2.2 times (FDR 2e-5) | 0.64 |
| Micro-PVM_4-SEAAD (disease-associated) | 26 | 18 | not significant | 1.32 |
| Micro-PVM_2 (homeostatic) | 62 | 215 | no, 7 times less common | 1.68 |

Disease-associated states are more common in AD donors, but they do not have higher risk scores
than homeostatic microglia.

### The Mic.12 / Mic.13 / Ast.10 test

| Test | Mic.12 | Mic.13 | Ast.10 |
|---|---|---|---|
| All nuclei, top 25 % vs rest | FDR 0.001 | FDR 2e-14 | FDR 0.91 |
| Inside own cell type, correlation (rho, p) | 0.18, 3e-6 | 0.14, 3e-4 | 0.02, 0.44 |
| Inside microglia, *APOE* removed, correlation (rho, p) | -0.03, 0.48 | 0.03, 0.37 | - |
| Inside microglia, *APOE* removed, top 25 % vs rest | p = 0.14 | p = 0.008 | - |

About 90 % of high-Mic.13 nuclei are not microglia, so the all-nuclei test mostly compares cell
types. Inside microglia the links are driven largely by *APOE*.

## Limitations

- **Small pilot.** 12 donors. AD status and donor cannot be separated, so differences in cell-state
  composition between AD and control may partly be donor differences.
- **Marker-gene proxies.** Mic.12, Mic.13 and Ast.10 are measured with 2 to 4 marker genes each,
  not with real state labels. *CPM*, *TREM2* and *MT1X* are missing from the dataset, and the
  microglial markers are also expressed by other cell types.
- **Missing genes.** The GEO matrix has only 10,850 genes. 75 of the 147 significant GWAS genes are
  not in it, including *TREM2*, *MS4A4E* and *ACE*. The *APOE* variants that define the e4 and e2
  alleles (rs429358, rs7412) are not in the GWAS file, so the *APOE* signal comes through nearby
  variants.
- **Different brain region.** This dataset is entorhinal cortex; Green et al. used prefrontal
  cortex.
- **No cell filtering.** The GEO matrix is already filtered by its authors; a standard QC filter
  would remove no nuclei, so none were removed.
- **scDRS group test reimplemented.** The group-level function in scdrs 1.0.2 does not run with
  current numpy/pandas. Notebook 07 reimplements its association test exactly; the heterogeneity
  statistic is not included.

## How to run

```bash
git clone https://github.com/ar0700288/AD-glial-states-genetic-risk.git
cd AD-glial-states-genetic-risk
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt
python run_pipeline.py
```

On Linux or macOS use `.venv/bin/pip`. Run single stages with `python run_pipeline.py 05 06 07`.
The full run takes about 30 to 35 minutes, mostly MAGMA (about 20 minutes).

Before running, download the input data ([data/README.md](data/README.md)) and MAGMA
([tools/README.md](tools/README.md)).

**Manual step.** The Allen MapMyCells labels come from a web tool and cannot be regenerated by code.
The result is included in `data/mapmycells/`. To redo it, upload `interim/mapmycells_raw_upload.h5ad`
(created by notebook 07) to MapMyCells and replace the CSV.

## Repository layout

```
data/raw/            downloaded input data (not in git)
data/mapmycells/     Allen MapMyCells result (in git)
notebooks/           01 to 08, run in order
src/adpipe.py        shared folder paths and save functions
interim/             large intermediate files, can be rebuilt (not in git)
results/tables/      result tables
results/figures/     result figures
tools/               MAGMA program and reference panel (not in git)
run_pipeline.py      runs all notebooks in order
```

## Result files

| File | Contents |
|---|---|
| `01_qc_summary.tsv` | QC summary for AD and control |
| `02_cluster_composition.tsv` | Leiden clusters vs author cell types |
| `03_gwas_qc_summary.tsv` | Variants removed by each GWAS filter |
| `04_magma_genes.tsv` | MAGMA risk Z score and p-value for 18,626 genes |
| `05_marker_coverage.tsv` | Which marker genes exist in the dataset |
| `07_cell_metadata.tsv` | One row per nucleus: all labels, risk scores and state scores |
| `07_mapmycells_concordance.tsv` | Author labels vs Allen labels |
| `07_scdrs_group_stats.tsv` | Risk test for every cell type, state and cluster |
| `07_microglia_supertypes.tsv` | Microglial states: AD vs control and risk scores |
| `08_state_enrichment.tsv` | Mic.12, Mic.13, Ast.10 risk tests |

Figures for each stage are in `results/figures/`.

## License

The code is released under the [MIT License](LICENSE). The datasets, MAGMA and the Allen reference
data keep their own terms of use.

## References

- Green GS, et al. Cellular communities reveal trajectories of brain ageing and Alzheimer's disease.
  *Nature* (2024). doi:10.1038/s41586-024-07871-6
- Grubman A, et al. A single-cell atlas of entorhinal cortex from individuals with Alzheimer's
  disease reveals cell-type-specific gene expression regulation. *Nat Neurosci* 22, 2087-2097 (2019).
- Bellenguez C, et al. New insights into the genetic etiology of Alzheimer's disease and related
  dementias. *Nat Genet* 54, 412-436 (2022).
- de Leeuw CA, et al. MAGMA: Generalized Gene-Set Analysis of GWAS Data. *PLoS Comput Biol* 11,
  e1004219 (2015).
- Zhang MJ, et al. Polygenic enrichment distinguishes disease associations of individual cells in
  single-cell RNA-seq data. *Nat Genet* 54, 1572-1580 (2022).
- Gabitto MI, et al. Integrated multimodal cell atlas of Alzheimer's disease. *Nat Neurosci* (2024).
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biol* 19, 15 (2018).
