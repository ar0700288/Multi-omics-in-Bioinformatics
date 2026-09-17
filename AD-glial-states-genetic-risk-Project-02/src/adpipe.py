"""Shared paths and helpers for the AD single-cell / GWAS pipeline.

Every notebook starts with::

    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path.cwd().parent / "src"))
    from adpipe import *

which gives it the project paths plus ``save_table`` / ``save_fig``.

The directory layout is deliberate:

* ``data/raw``      - downloaded inputs, never written to (git-ignored, large)
* ``data/mapmycells`` - the MapMyCells web-service result (committed; cannot be regenerated locally)
* ``interim``       - heavy intermediates that any notebook can rebuild (git-ignored)
* ``results/tables``, ``results/figures`` - the small, committed outputs

Only ``results/`` is meant to be read by a human or committed to git.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

__all__ = [
    "ROOT", "RAW", "MAPMYCELLS", "INTERIM", "MAGMA_DIR", "MAGMA_TOOL",
    "RESULTS", "TABLES", "FIGURES",
    "save_table", "save_fig", "run_magma", "describe_outputs",
]

# --------------------------------------------------------------------------- paths


def _find_root(start: Path | None = None) -> Path:
    """Walk up from `start` until the directory that holds `data/` and `notebooks/`."""
    here = (start or Path.cwd()).resolve()
    for cand in (here, *here.parents):
        if (cand / "data").is_dir() and (cand / "notebooks").is_dir():
            return cand
    raise RuntimeError(f"project root not found above {here}")


ROOT = _find_root()

RAW = ROOT / "data" / "raw"
MAPMYCELLS = ROOT / "data" / "mapmycells"
INTERIM = ROOT / "interim"
MAGMA_DIR = INTERIM / "magma"
MAGMA_TOOL = ROOT / "tools" / "magma"

RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

for _d in (INTERIM, MAGMA_DIR, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------- save helpers


def save_table(df: pd.DataFrame, name: str, index: bool = True) -> Path:
    """Write a small result table to ``results/tables`` and report what it contains."""
    if not name.endswith(".tsv"):
        name += ".tsv"
    path = TABLES / name
    df.to_csv(path, sep="\t", index=index)
    print(f"  table -> results/tables/{name}  ({len(df):,} rows x {df.shape[1]} cols, "
          f"{path.stat().st_size / 1024:.0f} KB)")
    return path


def save_fig(fig: plt.Figure, name: str, dpi: int = 150) -> Path:
    """Write a figure to ``results/figures``."""
    if not name.endswith(".png"):
        name += ".png"
    path = FIGURES / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"  figure -> results/figures/{name}  ({path.stat().st_size / 1024:.0f} KB)")
    return path


def describe_outputs() -> pd.DataFrame:
    """List everything currently in ``results/`` - used by the final notebook."""
    rows = []
    for p in sorted(RESULTS.rglob("*")):
        if p.is_file():
            rows.append({"file": str(p.relative_to(ROOT)).replace("\\", "/"),
                         "size_kb": round(p.stat().st_size / 1024, 1)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------- magma glue


def run_magma(args: list[str]) -> None:
    """Run the MAGMA binary from ``tools/magma`` (its reference panel is resolved
    relative to that directory) and stream its log."""
    exe = MAGMA_TOOL / ("magma.exe" if (MAGMA_TOOL / "magma.exe").exists() else "magma")
    if not exe.exists():
        raise FileNotFoundError(
            f"MAGMA binary not found at {exe}. See tools/README.md for how to install it."
        )
    cmd = [str(exe), *args]
    print("$", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=MAGMA_TOOL, capture_output=True, text=True)
    print(proc.stdout[-4000:])
    if proc.returncode != 0:
        print(proc.stderr[-4000:])
        raise RuntimeError(f"MAGMA failed with exit code {proc.returncode}")
