# NestBoot SCENIC+

Fifty paired bootstraps of enhancer-driven gene regulatory networks. Each pair reruns the stochastic GRNBoost2 TF-to-gene step and rebuilds eRegulons. The region-to-gene table is held fixed. The null shuffles expression independently within each gene before GRNBoost2.

The published numbers were produced with scenicplus 1.0a2 at commit `840dab85de3044846234c157e327b6ea0290abfb`. An earlier, unpinned eGRN pass wrote null networks one to two orders of magnitude smaller from the same TF-to-gene tables. Those files are `eRegulons.legacy.tsv` and are not used.

## Setup

Python 3.11.8 and [uv](https://astral.sh/uv). scenicplus 1.0a2 is not a PyPI release, and neither are the aertslab packages it was run with. Install those five git commits first, then the pins in `requirements.lock`. `setuptools` must stay below 81 because scenicplus still imports `pkg_resources`. dask 2024.2.1 breaks on CPython 3.11.9 and newer.

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

`requirements.txt` is only the scenicplus commit plus `mudata`, `pandas`, `numpy`, and `scipy`. That is enough for `nestboot_scenicplus.py stats` and the Figure 2 scripts. A GRNBoost2 rerun needs the lock above: NumPy 1.26.4, pandas 1.5.0, arboreto 0.1.6, dask and distributed 2024.2.1, and scikit-learn 1.3.2. On a CPU without AVX2, reinstall `polars-lts-cpu==0.20.13`. This machine did that.

## Inputs

| file | role |
|---|---|
| `work/nestboot_scenicplus/base_mudata.h5mu` | multiome, expression modality `scRNA` |
| `work/menr/tf_names.txt` | TF list for GRNBoost2 |
| `work/menr/cistromes_direct.h5ad` | direct cistromes for eGRN |
| `results_full_scenicplus/real_0/region_to_gene.tsv` | fixed region-to-gene table |
| `hg38_screen_v10_clust.regions_vs_motifs.rankings.feather` | cisTarget rankings used for motif enrichment |

`region_to_gene` is deterministic, but rebuilding it needs the SCENIC+ search-space table (`scenicplus grn_inference region_to_gene`). That search space is not in this bundle, so the benchmark starts from the table above.

## Run

From this directory, with the virtualenv active:

```bash
python nestboot_scenicplus.py run --n_cpu 40
python nestboot_scenicplus.py stats
```

`run` covers runs 0–49. A finished side (non-empty `tf_to_gene.tsv` and `eRegulons.tsv`) is skipped. A side that already has `tf_to_gene.tsv` repeats only eGRN. Seeds are `42 + run` for real and `42 + run + 10000` for the null. The null expression matrix is deleted after GRNBoost2.

`stats` rewrites `results_summary/` from the 50 `eRegulons.tsv` pairs. That step does not need the rankings database or GRNBoost2. It calls `build_nestboot_stats.py` and `nestboot_analysis.py`.

## Support, overlap, and frequency

These scripts read the finished `tf_to_gene.tsv` tables. They do not rerun GRNBoost2.

| script | result |
|---|---|
| `plot_nestboot.py` | eRegulon figures, including Equation 1 on TF–gene presence |
| `sparsity_fdr.py` | Equation 1 at 30 importance cutoffs |
| `raw_degree_fdr.py` | support tables for all edges and for 3–5 links per node, before dropping self-links |
| `inner_bootstrap.py` | 5×10 and 10×5 splits of the eRegulon links |
| `fig2_grid.py` | Figure 2 curves for 3–5 links per gene, 3–5 links per TF, and all edges, with self-links removed |
| `chapter_nestboot_fig2.py` | the same final figures, written to be read section by section |

`chapter_nestboot_fig2.py` writes `fig2_degree_3_and_5.png` and `fig2_all_edges.png`. Copies of those figures are in `results_summary/plots/`. The write-up is `results_summary/ANALYSIS.md`.

## What the statistics test

A TF–gene link or eRegulon name is present or absent in each run. Run `i` real is paired with run `i` null. The exact McNemar test asks whether disagreements fall on the real side. q-values are Benjamini–Hochberg within TF–gene links and within eRegulon names. See `results_summary/REPRODUCIBILITY.md` and `results_summary/NESTBOOT.md`.

## What is reproducible

The deposited eRegulon tables plus `python nestboot_scenicplus.py stats` reproduce the published counts and McNemar tests. A full rerun uses the deposited inputs and `python nestboot_scenicplus.py run`. Both depend on the environment in `requirements.lock` (Python 3.11.8, scenicplus commit `840dab85`, and the pinned NumPy, pandas, arboreto, dask, and scikit-learn versions). `scenicplus==1.0a2` on that list is not a PyPI release. Use the git installs in Setup.

This is as far as the released files go. They do not make a bitwise rerun:

- GRNBoost2 with a fixed seed can still change individual links at a different `--n_cpu`.
- The null shuffle uses NumPy's global random stream, which is part of the pinned NumPy 1.26.4 behavior.
- `region_to_gene.tsv` is an input. The SCENIC+ search space that built it is not in the deposit.
- The MuData, cistromes, and TF list are inputs. This repository does not rebuild them from raw multiome.
- Nineteen null eRegulon tables were rebuilt from saved `tf_to_gene.tsv` files after an earlier eGRN pass. The published files are the rebuilt ones.

## Zenodo

One unpublished draft, [23559287](https://zenodo.org/deposit/23559287). Inputs and the 50 real/null outputs are in that record. The run outputs are packed as `outputs_all.zip`. Filenames are flat because a Zenodo bucket key cannot contain directories. `zenodo_manifest.tsv` maps each key back to a relative path and records sha256.

```bash
python upload_zenodo.py --manifest-only
ZENODO_TOKEN=... python upload_zenodo.py --creator "Family, Given"
python restore_zenodo.py --src downloads/inputs --dest . --deposit inputs --check
python restore_zenodo.py --src downloads/outputs --dest . --deposit outputs --check
```

The token is a Zenodo personal access token with deposit write scope. The script does not publish the drafts. `--check` rereads every file to verify sha256.
