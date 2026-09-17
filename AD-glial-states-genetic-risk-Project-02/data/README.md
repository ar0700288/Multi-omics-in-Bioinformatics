# Data

## `raw/` — not tracked by git (≈770 MB)

Download these three files into `data/raw/` before running the pipeline.

### Single-nucleus RNA-seq — GEO [GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852)

Grubman A, Chew G, Ouyang JF, *et al.* "A single-cell atlas of entorhinal cortex from individuals
with Alzheimer's disease reveals cell-type-specific gene expression regulation."
*Nat Neurosci* 22, 2087–2097 (2019).

| File | Size |
|---|---|
| `GSE138852_counts.csv.gz` | 12 MB |
| `GSE138852_covariates.csv.gz` | 87 KB |

```bash
cd data/raw
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_counts.csv.gz
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_covariates.csv.gz
```

13,214 nuclei × 10,850 genes, 6 AD and 6 control entorhinal-cortex donors. The matrix is the
authors' post-QC set — the pipeline does not filter it further (see notebook 01).

### GWAS summary statistics — GWAS Catalog [GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158)

Bellenguez C, Küçükali F, Jansen IE, *et al.* "New insights into the genetic etiology of
Alzheimer's disease and related dementias." *Nat Genet* 54, 412–436 (2022).

| File | Size |
|---|---|
| `GCST90027158_buildGRCh38.tsv.gz` | 755 MB |

Download the harmonised GRCh38 file from the study page above. 111,326 clinically diagnosed /
proxy AD cases and 677,663 controls.

## `mapmycells/` — tracked by git (2 MB)

The result of mapping the nuclei against the Allen Institute **SEA-AD 10x Human MTG
(CCN20230505)** taxonomy with *Deep Generative Mapping*
(`vae_cell_type_mapper_scvi` v2.0.2).

This is committed because **it cannot be regenerated locally** — MapMyCells is a web service.
To reproduce it: run notebooks 01–02, then notebook 07's first cell writes
`interim/mapmycells_raw_upload.h5ad` (raw counts, cells × genes). Upload that at
<https://knowledge.brain-map.org/mapmycells/process>, choose the 10x Human MTG SEA-AD reference
with Deep Generative Mapping, and unpack the returned archive here, renaming the CSV to
`mapmycells_SEAAD_MTG_deepgenerative.csv`.

| File | What it is |
|---|---|
| `mapmycells_SEAAD_MTG_deepgenerative.csv` | class / subclass / supertype call + softmax confidence per nucleus |
| `summary_metadata.json` | 13,214 cells and 10,850 genes mapped |
| `validation_log.txt` | the service's input-validation report |
