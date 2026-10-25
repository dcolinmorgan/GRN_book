# Seed consensus for SCENIC+: chapter modules

The repository is [Scenic+/nestboot](https://github.com/dcolinmorgan/GRN_book/tree/main/Scenic+/nestboot). Four modules: install, load data, run, analysis. The runner is `nestboot_scenicplus.py`. The figures are `fig2_grid.py` and `chapter_nestboot_fig2.py`. Folder and file names keep the word "nestboot" from the original plan. The analysis below is a consensus across GRNBoost2 seeds with a shuffled control, and the last section explains why it is not NestBoot.

The Aerts-lab code is the five git commits in the install module. The cisTarget rankings database is a file in the Zenodo record, not a second download.

## Install

From `Scenic+/nestboot`, with Python 3.11.8 and uv:

```bash
uv venv --python 3.11.8 .venv
source .venv/bin/activate
uv pip install --no-deps \
  "scenicplus @ git+https://github.com/aertslab/scenicplus.git@840dab85de3044846234c157e327b6ea0290abfb" \
  "pycistarget @ git+https://github.com/aertslab/pycistarget@5aa517604e4842539a7531c16905825dc7cb80fb" \
  "pycisTopic @ git+https://github.com/aertslab/pycisTopic@787ce422a37f5975b0ebb9e7b19eeaed44847501" \
  "pyscenic @ git+https://github.com/aertslab/pySCENIC@06bafba412792f6efa5a552a23bb221cc3bdea1b" \
  "loomxpy @ git+https://github.com/aertslab/LoomXpy@61995ff10940968eac2cee8fe48300ab477a15d0"
grep -Ev '^(scenicplus|pycistarget|pycisTopic|pyscenic|loomxpy)==' requirements.lock > /tmp/requirements.pypi.txt
uv pip install --build-constraints build-constraints.txt -r /tmp/requirements.pypi.txt
uv pip install "setuptools<81"
```

`--no-deps` keeps `uv` from reading scenicplus's own requirement `pycistarget @ git+https://github.com/aertslab/pycistarget`, which has no commit and conflicts with the pin above. The same conflict exists for pycisTopic, pySCENIC, and LoomXpy. Drop those five names from the lock file before installing it. The git installs already provide them.

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

Each run changes only the GRNBoost2 seed. Cells are not resampled. It writes 50 pairs under `results_full_scenicplus/`. For run `i`:

1. Measured side, seed `42 + i`. GRNBoost2 (`scenicplus grn_inference TF_to_gene`, method GBM) scores TF–gene importance on the observed multiome, using 40 cores. The eGRN step then builds `eRegulons.tsv` from that table plus the fixed region-to-gene links, direct cistromes, and the rankings database.
2. Shuffled side, seed `42 + i + 10000`. Each gene’s expression is permuted across cells on its own. GRNBoost2 and the eGRN step are repeated on that matrix. The shuffled MuData file is deleted after GRNBoost2. The region-to-gene table, cistromes, and rankings are the same files as the measured side.

A side that already has both `tf_to_gene.tsv` and `eRegulons.tsv` is skipped. If only the TF–gene table is present, only the eGRN step is repeated. `--n_cpu 40` is the worker count for both scenicplus calls. The same seed at a different core count need not reproduce the deposited links.

## Analysis

```bash
python fig2_grid.py
```

That writes the support curves for 3, 4, and 5 links per gene, 3, 4, and 5 links per TF, and for every TF–gene edge, plus the table `results_summary/fig2_degree_and_full.tsv`. In each run, the degree cut keeps the strongest edges by GRNBoost2 importance. Self-links are removed first, because a TF predicting itself survives the shuffle. `chapter_nestboot_fig2.py` is the book version of the 3-and-5 figure and the all-edges figure.

### Seed support

GRNBoost2 is deterministic for a fixed matrix and a fixed seed. Changing the seed changes the network. A link’s seed support is the fraction of the 50 seeded runs that contain it. The curves show, for each support value, the fraction of links with exactly that support. They are not smoothed. The y-axis runs from 0 to 0.2, and the one-seed peak sits above that frame.

The measured networks recur far more than the shuffled ones. At 3 regulators per gene, 62% of measured links appear in at least 2 seeds, against 13% of shuffled links. Two measured runs share a Jaccard overlap of 0.30, against 0.003 for two shuffled runs. Among links in at least 40 of 50 seeds, the overlap rises to about 0.8. No shuffled link reaches 8 seeds. 924 measured links are in all 50.

### Cutoff from the shuffled control

The cutoff is set by the shuffled runs, not by a chosen percentage. For each support B, count the shuffled links and the measured links with support above B. The cutoff is the smallest B at which

shuffled links above B / measured links above B ≤ 0.05.

That ratio is an empirical false-discovery rate against the shuffled control. The black line in each panel marks it.

| Graph | Kept when seen in at least | Measured links | Shuffled links | Ratio |
|---|---:|---:|---:|---:|
| 3 TFs per gene | 4 of 50 seeds | 230,264 | 4,069 | 0.018 |
| 5 TFs per gene | 5 of 50 seeds | 320,337 | 3,391 | 0.011 |
| 3 targets per TF | 2 of 50 seeds | 13,775 | 683 | 0.050 |
| 5 targets per TF | 3 of 50 seeds | 15,111 | 19 | 0.001 |
| All edges | 22 of 50 seeds | 1,566,575 | 55,314 | 0.035 |

The cutoff is low on the degree-cut graphs because a shuffled matrix almost never produces the same edge twice. The set it keeps is therefore a union across seeds, larger than any single run: 230,264 links against 87,495 edges per run at 3 regulators per gene. Report seed support with each link, so a reader can tell the links that pass the control from the core that recurs in most seeds. With every edge kept, shuffled links persist up to support 0.78, and the cutoff moves out to 22 of 50 seeds. The degree cut removes that tail.

### Why this is not NestBoot

NestBoot (Morgan et al., *Bioinformatics* 35:1026–1032, 2019) bootstraps independent experiments in an inner loop, repeats that inner loop in an outer loop with a fresh shuffle, and keeps a link only when it clears the Equation 1 cutoff in every outer run. It assumes the inference is fixed once the input is fixed, so all variation comes from the data. The ratio used above is the paper’s Equation 1.

None of that applies to these runs. There is one multiome, so there are no experiments to resample. The measured matrix is the same in all 50 runs, and only the GRNBoost2 seed changes. The variation is the solver’s, not the sample’s, so there is one set of runs and no outer loop. Splitting the 50 seeds into groups, such as 5 groups of 10, only repartitions the same seeds. The shuffled side permutes each gene within the one dataset rather than shuffling experiments.

A real NestBoot here would draw cells with replacement, keep the GRNBoost2 seed fixed, and repeat Q inner and R outer times on measured and shuffled data. One GRNBoost2 fit took about 1.5 hours on this dataset, so even R = 5 and Q = 20 is 200 fits. That was not run.
