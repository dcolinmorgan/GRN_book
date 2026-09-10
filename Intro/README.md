# Introduction Chapter — Figures

Reproducible code for the figures in the introductory chapter,
*"Biological Insights from Gene Regulatory Networks: An Overview of
Approaches for Inference and Analysis."*

## Requirements

```bash
pip install numpy pandas matplotlib scikit-learn
```

## Figure 2 — Dataset landscape (`plot_datasets.py`)

A scatter plot surveying datasets used for GRN inference, positioned by genes
perturbed (x) vs. total cells/samples profiled (y, log scale). Encodings:

- **Marker size** — genes measured (median genes/cell for droplet single-cell
  studies; annotated feature space for bulk)
- **Marker shape** — organism
- **Fill color** — perturbation type (concentric rings for mixed-type datasets
  such as L1000)

```bash
python plot_datasets.py        # writes datasets_scatter.png / .pdf
```

The dataset table is curated inline in the script (from a source TSV plus
primary-source corrections); edit the `datasets` DataFrame to add or revise
entries.

## Figure 4 — Why benchmarking metrics mislead under class imbalance (`plot_grn_reliability.py`)

A five-panel figure built on a **toy model with no explicit topology**: over the
set of possible directed edges, each edge is labeled true/false (mean out-degree
~6, i.e. `6·N` true edges) and given a hypothetical inferred score drawn from two
overlapping Gaussians — true edges from `N(1.5, 1)`, non-edges from `N(0, 1)`.
The distribution separation (d = 1.5σ) is held constant across all network
sizes, so scorer *quality* is fixed; only class prevalence changes.

- **A** — prevalence (= true/possible edges ≈ k/N) falls as networks grow
- **B** — ROC curves are nearly identical across sizes (AUROC is
  prevalence-invariant)
- **C** — PR curves collapse toward each size's prevalence baseline (dashed)
- **D** — *same scorer, different verdicts:* AUROC stays flat, AUPR collapses,
  and two normalization attempts are contrasted (see below)
- **E** — negative subsampling inflates a poor AUPR at every network size

```bash
python plot_grn_reliability.py   # writes grn_reliability_metrics.png / .pdf + per-panel TSVs
```

### On the DREAM5-style permutation z-score (Panel D)

Raw AUPR shrinks with network size (driven by falling prevalence), and the naive
**AUPR-ratio** (AUPR / prevalence) *over-corrects* and rises with size, because a
genuinely discriminative scorer's AUPR decays sub-linearly in prevalence. Neither
is size-invariant.

The more robust option, implemented in `permutation_z()`, is a **DREAM5-style
permutation z-score**: build an empirical null by shuffling the true/false labels
many times and recomputing AUPR, then report

```
z = (observed_AUPR − mean(null_AUPR)) / std(null_AUPR)
```

To make the z-score comparable *across* network sizes, every size is evaluated on
a **fully matched subsample** (identical positive and negative counts), so both
the evaluation prevalence and the null resolution are the same. This keeps the
z-score roughly flat for a fixed-quality scorer.

Caveats (see also Badia-i-Mompel et al., *Nat Rev Genet* 2023, and DREAM
challenge literature):

- No AUPR-derived statistic is *provably* scale-invariant; the matched-null
  z-score is an empirical fix, not an algebraic one.
- Networks too small to supply the matched subsample (here roughly N < 50)
  skew the z-score downward — their possible-negative pool is too small to
  resolve the null.
- More permutations (300 → 2000) barely move the estimate; the null mean/std
  converge quickly, so instability is a matched-size issue, not a
  too-few-randomizations issue.

Per-panel data are also written as `grn_reliability_panel*.tsv` for independent
re-plotting.
