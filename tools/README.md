# External tools

Notebook 04 shells out to MAGMA. The binary and its LD reference panel are not tracked by git
(≈600 MB, and MAGMA has its own licence). Install them into `tools/magma/` before running
stage 04.

## MAGMA v1.10

Download from <https://cncr.nl/research/magma/> — pick the build for your platform
(`magma_v1.10_win.zip`, `magma_v1.10.zip` for Linux, `magma_v1.10_mac.zip`) and unzip into
`tools/magma/`.

## Gene locations — NCBI38

From the same page, download `NCBI38.zip` and unzip `NCBI38.gene.loc` into `tools/magma/`.

## LD reference panel — 1000 Genomes European

From the same page, download `g1000_eur.zip` and unzip into `tools/magma/`. This gives
`g1000_eur.bed`, `.bim`, `.fam` and `.synonyms` (503 individuals, 22.6 M variants).

## Expected layout

```
tools/magma/
├── magma.exe            (or `magma` on Linux/macOS)
├── NCBI38.gene.loc
├── g1000_eur.bed
├── g1000_eur.bim
├── g1000_eur.fam
└── g1000_eur.synonyms
```

`src/adpipe.py::run_magma` runs the binary with `tools/magma/` as its working directory, which is
why `--bfile g1000_eur` and `--gene-loc NCBI38.gene.loc` are given as bare names.

## Citation

de Leeuw CA, Mooij JM, Heskes T, Posthuma D. "MAGMA: Generalized Gene-Set Analysis of GWAS Data."
*PLoS Comput Biol* 11(4): e1004219 (2015).
