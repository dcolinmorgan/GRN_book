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
  and the AUPR-ratio normalization is contrasted (see below)
- **E** — negative subsampling inflates a poor AUPR at every network size

```bash
python plot_grn_reliability.py   # writes grn_reliability_metrics.png / .pdf + per-panel TSVs
```

### On AUPR normalization and the DREAM5 significance transform (Panel D)

Raw AUPR shrinks with network size (driven by falling prevalence). The common
**AUPR-ratio** (AUPR / prevalence, i.e. AUPR over the expected-random baseline)
is meant to divide out that floor, but it *over-corrects* and RISES with size,
because a genuinely discriminative scorer's AUPR decays sub-linearly in
prevalence. So it is not size-invariant either.

A natural next idea is a **DREAM5-style significance transform** (Marbach et
al., *Nat Methods* 2012): build an empirical null of AUPR by shuffling the
true/false labels, then report a permutation z-score or `-log10(p)` of the
observed AUPR against that null. DREAM5 used p-values (with a stretched-
exponential tail fit), not z-scores. The catch — confirmed numerically here — is
that this is a **significance** measure, not an effect size: with thousands to
millions of edges the null becomes extremely tight, so any non-random scorer
yields astronomically small p-values (z-scores of 30+, p far below 1e-300).
These numbers are dominated by sample size, not prediction quality, and are
therefore not a bounded, size-invariant quality metric.

**Conclusion:** no AUPR-derived scalar (ratio, permutation z-score, or
`-log10(p)`) gives a bounded, size-invariant quality number. **AUROC** is the
only prevalence-invariant summary shown, though it is comparatively insensitive
to imbalance (Panel B). Cross-size comparison of AUPR is only meaningful at
matched prevalence / negative-sampling ratio (Panel E).

Per-panel data are also written as `grn_reliability_panel*.tsv` for independent
re-plotting.
