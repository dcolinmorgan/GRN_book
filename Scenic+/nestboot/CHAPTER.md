# NestBoot SCENIC+ chapter modules

The repository is [Scenic+/nestboot](https://github.com/dcolinmorgan/GRN_book/tree/main/Scenic+/nestboot). Four modules: install, load data, run, analysis. The runner is `nestboot_scenicplus.py`. The NestBoot figures are `fig2_grid.py` and `chapter_nestboot_fig2.py`.

The Aerts-lab code is the five git commits in the install module. The cisTarget rankings database is a file in the Zenodo record, not a second download.

## Install

From `Scenic+/nestboot`, with Python 3.11.8 and uv:

```bash
uv venv --python 3.11.8 .venv
source .venv/bin/activate
uv pip install --build-constraints build-constraints.txt \
  "scenicplus @ git+https://github.com/aertslab/scenicplus.git@840dab85de3044846234c157e327b6ea0290abfb" \
  "pycistarget @ git+https://github.com/aertslab/pycistarget@5aa517604e4842539a7531c16905825dc7cb80fb" \
  "pycisTopic @ git+https://github.com/aertslab/pycisTopic@787ce422a37f5975b0ebb9e7b19eeaed44847501" \
  "pyscenic @ git+https://github.com/aertslab/pySCENIC@06bafba412792f6efa5a552a23bb221cc3bdea1b" \
  "loomxpy @ git+https://github.com/aertslab/LoomXpy@61995ff10940968eac2cee8fe48300ab477a15d0"
uv pip install --build-constraints build-constraints.txt -r requirements.lock
uv pip install "setuptools<81"
```

`requirements.txt` is enough for the statistics and the figures. The lock is required for a GRNBoost2 rerun. On a CPU without AVX2, install `polars-lts-cpu==0.20.13` instead of `polars`.

## Load data

One public Zenodo record holds the inputs and the 50 outputs: [10.5281/zenodo.23559287](https://doi.org/10.5281/zenodo.23559287). The rankings database is in that record.

```bash
zenodo_get 10.5281/zenodo.23559287
```

`zenodo_get` reads the public record API and saves every file. Unpack the runs with `unzip outputs_all.zip`. The downloaded names are flat (`work__nestboot_scenicplus__base_mudata.h5mu`, and so on). `restore_zenodo.py` belongs to an older two-deposit layout. Do not use it for this record.

To recompute the figures from the deposited runs, stop here and go to Analysis. The rankings file is only required for a fresh eGRN inference.

## Run

```bash
python nestboot_scenicplus.py run --n_cpu 40
```

This does not resample cells and it does not repeat an outer loop. It writes 50 pairs under `results_full_scenicplus/`. For run `i`:

1. Measured side, seed `42 + i`. GRNBoost2 (`scenicplus grn_inference TF_to_gene`, method GBM) scores TF–gene importance on the observed multiome, using 40 cores. The eGRN step then builds `eRegulons.tsv` from that table plus the fixed region-to-gene links, direct cistromes, and the rankings database.
2. Shuffled side, seed `42 + i + 10000`. Each gene’s expression is permuted across cells on its own. GRNBoost2 and the eGRN step are repeated on that matrix. The shuffled MuData file is deleted after GRNBoost2. The region-to-gene table, cistromes, and rankings are the same files as the measured side.

A side that already has both `tf_to_gene.tsv` and `eRegulons.tsv` is skipped. If only the TF–gene table is present, only the eGRN step is repeated. `--n_cpu 40` is the worker count for both scenicplus calls. The same seed at a different core count need not reproduce the deposited links.

## Analysis

NestBoot, at one sparsity, scores each link by how often the 50 seeds infer it. Support is that fraction. Equation 1 sets the cutoff B:

FDR(B) = (shuffled links with support > B) / (measured links with support > B)

The cutoff is the smallest B at which this ratio is at most 0.05. The curves are not smoothed. Each point is the fraction of links at that exact support.

```bash
python fig2_grid.py
```

That writes the frequency curves for 3, 4, and 5 links per gene, 3, 4, and 5 links per TF, and for every TF–gene edge. Self-links are removed before the cut. `chapter_nestboot_fig2.py` is the book version of the 3-and-5 figure and the all-edges figure. The y-axis is frequency from 0 to 0.2. The one-run peak sits above that frame.

These 50 seeds are one inner loop. They are not 100 outer runs of experiment bootstraps, so a link is kept when its support exceeds B. There is no further “pass every outer run” filter.
