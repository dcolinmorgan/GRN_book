# NestBoot SCENIC+ chapter modules

The repository is [Scenic+/nestboot](https://github.com/dcolinmorgan/GRN_book/tree/main/Scenic+/nestboot). Four modules: install, load data, run, analysis. The runner is `nestboot_scenicplus.py`. `nestboot_analysis.py` is not a separate entry point. `python nestboot_scenicplus.py stats` calls it.

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

Runs 0–49. A finished side is skipped. Seeds are `42 + run` for the measured data and `42 + run + 10000` for the shuffle. The shuffle permutes each gene across cells. Region-to-gene, cistromes, and the rankings database stay fixed. GRNBoost2 at a different `--n_cpu` need not match the deposited `tf_to_gene.tsv` link for link.

## Analysis

eRegulon counts and the McNemar tests:

```bash
python nestboot_scenicplus.py stats
```

That calls `build_nestboot_stats.py` and `nestboot_analysis.py`. It reads `eRegulons.tsv` only.

Link frequency at 3 and 5 links per node, and for every TF–gene edge, with self-links removed:

```bash
python fig2_grid.py
```

`chapter_nestboot_fig2.py` is the same frequency figure written section by section for the book text. The curves are not smoothed. Each point is the fraction of links at that exact support.

Equation 1 at FDR 0.05 is the smallest support B where the number of shuffled links above B, divided by the number of measured links above B, is at most 0.05. These 50 seeds are one inner loop. They are not 100 outer runs of experiment bootstraps.
