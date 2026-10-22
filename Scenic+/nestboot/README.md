# NestBoot SCENIC+

Fifty paired bootstraps of enhancer-driven gene regulatory networks. Each pair reruns the stochastic GRNBoost2 TF-to-gene step and rebuilds eRegulons. The region-to-gene table is held fixed. The null shuffles expression independently within each gene before GRNBoost2.

The published numbers were produced with scenicplus 1.0a2 at commit `840dab85de3044846234c157e327b6ea0290abfb`. An earlier, unpinned eGRN pass wrote null networks one to two orders of magnitude smaller from the same TF-to-gene tables. Those files are `eRegulons.legacy.tsv` and are not used.

## Setup

Python 3.11.8 and [uv](https://astral.sh/uv).

```bash
uv venv --python 3.11.8 .venv
source .venv/bin/activate
uv pip install --build-constraints build-constraints.txt -r requirements.txt
uv pip install "setuptools<81"
```

On a CPU without AVX2, reinstall the matching `polars-lts-cpu` wheel. `setup_new_machine.sh` does that check.

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
python nestboot.py run --n_cpu 40
python nestboot.py stats
```

`run` covers runs 0–49. A finished side (non-empty `tf_to_gene.tsv` and `eRegulons.tsv`) is skipped. A side that already has `tf_to_gene.tsv` repeats only eGRN. Seeds are `42 + run` for real and `42 + run + 10000` for the null. The null expression matrix is deleted after GRNBoost2.

`stats` rewrites `results_summary/` from the 50 `eRegulons.tsv` pairs. That step does not need the rankings database or GRNBoost2.

## What the statistics test

A TF–gene link or eRegulon name is present or absent in each run. Run `i` real is paired with run `i` null. The exact McNemar test asks whether disagreements fall on the real side. q-values are Benjamini–Hochberg within TF–gene links and within eRegulon names. See `results_summary/REPRODUCIBILITY.md` and `results_summary/NESTBOOT.md`.

## What is reproducible

The deposited eRegulon tables plus `python nestboot.py stats` reproduce the published counts and McNemar tests. A full rerun uses the deposited inputs and `python nestboot.py run`. Both depend on the environment in `requirements.lock` (Python 3.11.8, scenicplus commit `840dab85`, and the pinned NumPy, pandas, arboreto, dask, and scikit-learn versions). Install scenicplus from that git commit first; `scenicplus==1.0a2` on the lock list is not a PyPI release.

This is as far as the released files go. They do not make a bitwise rerun:

- GRNBoost2 with a fixed seed can still change individual links at a different `--n_cpu`.
- The null shuffle uses NumPy's global random stream, which is part of the pinned NumPy 1.26.4 behavior.
- `region_to_gene.tsv` is an input. The SCENIC+ search space that built it is not in the deposit.
- The MuData, cistromes, and TF list are inputs. This repository does not rebuild them from raw multiome.
- Nineteen null eRegulon tables were rebuilt from saved `tf_to_gene.tsv` files after an earlier eGRN pass. The published files are the rebuilt ones.

## Zenodo

Two drafts. Filenames are flat because a Zenodo bucket key cannot contain directories. `zenodo_manifest.tsv` maps each key back to a relative path and records sha256.

- inputs: MuData, TF names, cistromes, region-to-gene, rankings database, and the runner
- outputs: 50 real and 50 null `eRegulons.tsv` and `tf_to_gene.tsv`, `results_summary/`, and the runner

```bash
python upload_zenodo.py --manifest-only
ZENODO_TOKEN=... python upload_zenodo.py --creator "Family, Given"
python restore_zenodo.py --src downloads/inputs --dest . --deposit inputs --check
python restore_zenodo.py --src downloads/outputs --dest . --deposit outputs --check
```

The token is a Zenodo personal access token with deposit write scope. The script does not publish the drafts. `--check` rereads every file to verify sha256.
