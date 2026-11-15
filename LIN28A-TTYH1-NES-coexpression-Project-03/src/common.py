"""Settings and helpers shared by all steps.

Every notebook in notebooks/ imports from here, so paths, gene names, age bins and plot style are defined once.
All values come from config.yaml; change the config rather than the scripts.
"""
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# The notebooks first move to the repository root, so these relative paths work.
CFG = yaml.safe_load(open("config.yaml"))
B = CFG["braun"]                       # dataset settings: file path and metadata column names
GENES = CFG["genes_of_interest"]       # LIN28A, TTYH1, NES
BINS = CFG["age_bins"]["labels"]       # early, mid, late

# Intermediate files go to data/processed (not committed), final figures and tables to results/
PROCESSED = Path("data/processed")
TABLES = Path("results/tables")
FIGURES = Path("results/figures")
for folder in (PROCESSED, TABLES, FIGURES):
    folder.mkdir(parents=True, exist_ok=True)

# Same two colours for the 10x kits in every figure, so v2 and v3 are easy to follow
KIT_COLORS = {"v2": "#2a78d6", "v3": "#eb6834"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.dpi": 200, "savefig.bbox": "tight"})


def open_braun():
    # Backed mode: only the cell and gene tables are loaded into memory. The 1.67 million x 59,459
    # count matrix stays on disk and is read in blocks of cells by the scripts that need it.
    return ad.read_h5ad(B["file"], backed="r")


def cell_table(adata):
    """One row per cell with the metadata used in the analysis.

    Age and total UMIs are stored as text categories in the file, so they are converted to numbers.
    Ages are rounded to one decimal so that, for example, 6.6 and 6.60001 fall in the same week.
    """
    o = adata.obs
    cells = pd.DataFrame({
        "cell_class": o[B["class_col"]].astype(str),   # authors' annotation, e.g. Radial glia
        "cluster": o["cluster_id"].astype(str),        # authors' fine clusters (used to define the population)
        "region": o["Region"].astype(str),
        "chemistry": o[B["chemistry_col"]].astype(str),  # 10x kit: v2 or v3
        "donor": o[B["donor_col"]].astype(str),        # one embryo or fetus
        "age": pd.to_numeric(o[B["age_col"]].astype(str)).round(1),  # weeks post-conception
        "umis": pd.to_numeric(o["total_UMIs"].astype(str)),          # sequencing depth of the cell
    })
    # Age bins from config: early 5-7, mid 7.5-9.5, late 10-14 weeks
    cells["age_bin"] = pd.cut(cells["age"], CFG["age_bins"]["edges"], labels=BINS).astype(str)
    return cells


def save(fig, name):
    # Save the figure as PNG into results/figures, then show it in the notebook (show also closes it)
    fig.savefig(FIGURES / f"{name}.png")
    plt.show()
