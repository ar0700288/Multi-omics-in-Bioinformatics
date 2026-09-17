#!/usr/bin/env python
"""Execute the notebook pipeline end to end, in order, writing outputs back in place.

    python run_pipeline.py              # run every notebook
    python run_pipeline.py 05 06 07     # run only these stages
    python run_pipeline.py --list       # show the stages

Each notebook is executed with its working directory set to ``notebooks/`` so the relative
paths inside behave exactly as they do interactively. Execution stops at the first cell that
raises, and the traceback is printed.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parent
NOTEBOOK_DIR = ROOT / "notebooks"

# Per-stage timeout in seconds. MAGMA's gene analysis is the long one.
TIMEOUTS = {"04": 5400}
DEFAULT_TIMEOUT = 3600


def stages() -> list[Path]:
    return sorted(NOTEBOOK_DIR.glob("[0-9][0-9]_*.ipynb"))


def run_one(path: Path) -> tuple[bool, float]:
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=TIMEOUTS.get(path.name[:2], DEFAULT_TIMEOUT),
        kernel_name="python3",
        resources={"metadata": {"path": str(NOTEBOOK_DIR)}},
        allow_errors=False,
    )
    t0 = time.time()
    try:
        client.execute()
        ok = True
    except CellExecutionError as exc:
        print(f"\n!! {path.name} failed:\n{exc}", file=sys.stderr)
        ok = False
    finally:
        nbformat.write(nb, path)
    return ok, time.time() - t0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("only", nargs="*", help="stage number prefixes to run, e.g. 05 06")
    ap.add_argument("--list", action="store_true", help="list the stages and exit")
    args = ap.parse_args()

    todo = stages()
    if args.list:
        for p in todo:
            print(" ", p.name)
        return 0
    if args.only:
        todo = [p for p in todo if p.name[:2] in set(args.only)]
        if not todo:
            print("no stage matched", args.only, file=sys.stderr)
            return 2

    print(f"running {len(todo)} stage(s) from {NOTEBOOK_DIR}\n")
    total = 0.0
    for path in todo:
        print(f"--> {path.name}", flush=True)
        ok, secs = run_one(path)
        total += secs
        print(f"    {'ok' if ok else 'FAILED'} in {secs:.1f}s\n", flush=True)
        if not ok:
            return 1
    print(f"pipeline complete in {total / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
