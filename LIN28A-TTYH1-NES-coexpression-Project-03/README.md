# LIN28A, TTYH1 and NES co-expression in the developing human brain

Take-home assignment, Gargiulo Lab (MDC Berlin).

1. In which cell types or states of human embryonic and fetal development are `LIN28A`, `TTYH1` and
   `NES` co-expressed, judged cell by cell?
2. Which genes distinguish that population from other cell types at the same stage (cell-type axis),
   from the same lineage at a later stage (stage axis), or both?

## Data

| | |
|---|---|
| Dataset | Braun et al. 2023, *Comprehensive cell atlas of the first-trimester developing human brain*, Science 382, eadf1226, [doi:10.1126/science.adf1226](https://doi.org/10.1126/science.adf1226) |
| File | `human_dev.h5ad` from [github.com/linnarsson-lab/developing-human-brain](https://github.com/linnarsson-lab/developing-human-brain); download URL in `config.yaml` |
| Raw reads | EGA, EGAS00001004107 (controlled access, not used) |
| Content | 1,665,937 cells, 59,459 genes, raw UMI counts, 26 donors, 5 to 14 weeks post-conception, 10x v2 and v3 |

## How to run

```bash
conda env create -f environment.yml
conda activate bif01
snakemake --cores 1
```

Each step is a Jupyter notebook in `notebooks/`, split into sections. Snakemake runs every notebook
from top to bottom and saves it in place, so after a run each notebook also shows its printed output and
figures. The notebooks in this repository already include detailed explanations of every cell and the
outputs of the run of 13 November 2026. To work step by step instead, open the notebooks in Jupyter or VS Code and run `01` to `07` in
order. The count matrix is never loaded whole. It is read in blocks of `chunk_size` cells (`config.yaml`), so
the pipeline runs on a 16 GB laptop.

| Step | What it does | Output |
|---|---|---|
| `01_download.ipynb` | Downloads the 12 GB h5ad | `data/raw/` |
| `02_qc.ipynb` | Missing values in counts and metadata; batch effect of the 10x kit | `02_qc_batch.png` |
| `03_data_overview.ipynb` | Confirms raw counts and the three genes; cells, donors and kit per week | `03_data_overview.png` |
| `04_gene_detection.ipynb` | Reads the three genes for all cells; detection per class, age and kit | `04_gene_detection.png`, gene counts |
| `05_coexpression.ipynb` | Q1: co-expression per cell, per donor; selects the population | three `05_*.png`, population clusters |
| `06_pseudobulk.ipynb` | Sums counts per state, donor and kit | `pseudobulk.npz` |
| `07_markers.ipynb` | Q2: DESeq2 on both axes | three ranked tables, two `07_*.png` |

Figures go to `results/figures/`, marker tables to `results/tables/`. Intermediate files go to
`data/processed/` and are not committed.

The analysis has no random steps (no subsampling, no random initialisation), so no seed is needed and
a rerun gives the same numbers. Package versions are pinned in `environment.yml`.

## Repository layout

```
config.yaml            all settings: data link, genes, age bins, thresholds
Snakefile              runs the notebooks 01 to 07 in order
environment.yml        conda environment
notebooks/01_ to 07_   one notebook per step, one cell per section
src/common.py          shared paths, helpers and plot style, imported by every notebook
results/figures/       PNG figures
results/tables/        ranked marker tables
```

## Methods

**Data and quality control.** The authors' quality-controlled cells and cell-class labels are used
as published. `X` holds raw integer UMI counts with no missing or negative values. The metadata used
(cell class, cluster, region, kit, donor, age, total UMIs) has no gaps. Most weeks have one or two
donors, so ages are grouped: early (5 to 7 weeks), mid (7.5 to 9.5), late (10 to 14).

**Batch effect.** Radial glia counts summed per donor and kit were compared by PCA. PC1 (31% of
variation) follows age. PC2 (25%) separates the two 10x kits completely. v3 cells also have about 2.5
times more UMIs than v2 cells. Age is the biology of interest, the kit is a technical batch and the two
overlap (mid mostly v2, late mostly v3). The counts are not corrected. The kit is handled in the
analysis instead: detection and co-expression are shown per kit, the cell-type test is paired within
donor (kit belongs to the donor, so it cancels out) and the stage test has the kit as a covariate.

To check whether one kit term is enough, the kit effect per gene (v3 minus v2) was measured in the
early bin and compared with the mid and late bins. The early bin has 6 v2 and 6 v3 donors of the same age.
The correlation is 0.72 for mid and 0.53 for late. After removing the early kit effect, PC1 (37%)
follows age. Early v2 and v3 samples mix, which is expected since the correction was measured on them.
The informative bins are the others: the mid v3 donor moves in among the mid v2 donors of similar age,
while late v3 donors still sit apart from the single late v2 donor. So in the late bin kit and age
cannot be fully separated.

**Expression and co-expression.** A gene is expressed in a cell if it has at least 1 UMI. A cell
co-expresses the three genes (triple-positive) if all three are expressed. For every cell class x age
bin x kit x donor, the number of triple-positive cells is compared with the number expected if the
genes were independent, n x p(LIN28A) x p(TTYH1) x p(NES). The expectation is computed within 5 depth
groups of cells, so deeply sequenced cells cannot create co-expression by detecting more of
everything. Donors are the replicates.

**Population.** LIN28A has at most 10 UMIs in any cell, so whether one cell reads triple-positive is
largely chance. The population is therefore defined on the authors' clusters: clusters with at least
500 v3 cells and a triple-positive share at least 5 times the dataset-wide v3 share (0.235%).

**Markers.** Counts are summed per state x donor x kit (pseudobulk) and tested with DESeq2
(PyDESeq2). Cell-type axis: population vs all other early cell classes except radial glia (the population is
itself radial glia), paired by donor
(`~pair + group`). Stage axis: population vs late radial glia, the same lineage, all donors with the
kit as covariate (`~chemistry + group`). A gene is up in the population with padj < 0.05 and log2 fold change >= 1. Tables are
ranked by the Wald statistic and give log2 fold change, its standard error, p, padj and the % of cells
expressing the gene in each group.

## Results

**Co-expression.** 2,395 cells are triple-positive, 81% of them radial glia. TTYH1 and NES are
detected in most radial glia (94% and 76%, early v3). LIN28A is detected in about 1%, so it decides
which cells co-express all three. Triple-positive cells are most frequent in early neural crest,
radial glia (median 0.59% per donor, 6 donors) and intermediate progenitors. They fall to about 0.1% in
late radial glia.

**Population.** 8 radial glia clusters from telencephalon, diencephalon and early whole-brain
samples, median age 6 to 6.9 weeks, 50,799 cells.

**Markers.**

| Comparison | Samples | Up in population |
|---|---|---|
| Cell-type axis | 11 vs 11 | 2,443 |
| Stage axis (kit in model) | 11 vs 5 | 1,796 |
| Both axes | | 455 |

- Cell type only: radial glia genes that stay on later, such as TTYH1, SHROOM3 and VIM.
- Stage only: LIN28B, LIN28B-AS1, TRIM71, BCAT1 and LIN28A itself.
- Both: SALL4, DPPA4, TERT, HMGA1, LIX1, CCND1, NES, LRP2, LGR5, HES4, HES5, WNT4, WNT7B, FZD10.

## Limitations

- The dataset starts at week 5. LIN28A is already declining, so its peak is probably earlier than
  the data covers.
- LIN28A is sparse (at most 10 UMIs per cell), so single-cell co-expression rates are lower bounds.
- No donor was sequenced with both kits, so the kit effect is estimated between donors and includes
  some donor-to-donor variation.
- One or two donors per week. The late bin has only one v2 donor, so the kit effect in the stage test
  is estimated mostly from the early donors. The kit effect is also only partly the same at late ages
  (r = 0.53). Stage markers with small fold changes should be read with this in mind.
- Early samples are often whole-brain dissections and late ones regional, so some stage markers may
  partly reflect region.
- Cell classes and clusters come from the authors; the population thresholds are a choice set in
  `config.yaml`.

## Use of AI tools

An AI assistant (Claude) helped write and comment the code and draft text. I made the analytical
decisions, checked every result and can explain each step.
