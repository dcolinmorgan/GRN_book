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

A six-panel figure built on a **toy model with no explicit topology**: over the
set of possible directed edges, each edge is labeled true/false (mean out-degree
~6, i.e. `6·N` true edges) and given a hypothetical inferred score drawn from two
overlapping Gaussians — true edges from `N(1.5, 1)`, non-edges from `N(0, 1)`.
The distribution separation (d = 1.5σ) is held constant across all network
sizes, so scorer *quality* is fixed; only class prevalence changes.

- **A** — prevalence (= true/possible edges ≈ k/N) falls as networks grow
- **B** — ROC curves are nearly identical across sizes (AUROC is
  prevalence-invariant)
- **C** — PR curves collapse toward each size's prevalence baseline (dashed)
- **D** — *same scorer, different verdicts:* AUROC (green) stays flat, while both
  AUPR (red) and MCC at the top-k operating point (purple) decline with size —
  showing that even MCC is not fully imbalance-proof; only AUROC is invariant
- **E** — negative subsampling inflates a poor AUPR at every network size
  (log-scale sweep over neg:pos ratio)
- **F** — the plain contrast: reported AUPR on a **balanced 1:1 subsample** vs.
  on the **full edge set**, same scorer. Balanced subsampling reports ~0.85 at
  every size while the honest full-set AUPR collapses (0.54 → 0.13 → 0.02)

```bash
python plot_grn_reliability.py   # writes grn_reliability_metrics.png / .pdf + per-panel TSVs
```

### On metric choice under class imbalance

Raw AUPR shrinks with network size (driven by falling prevalence). The common
**AUPR-ratio** (AUPR / prevalence) is meant to divide out the random floor, but
it *over-corrects* and rises with size, because a discriminative scorer's AUPR
decays sub-linearly in prevalence — so it is not size-invariant either.

**MCC** (Matthews correlation coefficient), shown in Panel D at the top-k
operating point, folds in true negatives and is more balanced than AUPR, but in
this setting it still declines with size — it is more robust, not invariant.

A **DREAM5-style significance transform** (Marbach et al., *Nat Methods* 2012)
is sometimes suggested: build an empirical null of AUPR by shuffling labels and
report a p-value (DREAM5 used p-values with a stretched-exponential tail fit,
combined only to pool across sub-challenges — not z-scores, and not a cross-size
quality score). We verified numerically that with realistic edge counts the null
is so tight that any non-random scorer yields astronomically small p-values
(`-log10(p)` saturates), so it is a **significance** measure, not a bounded
quality metric.

**Conclusion:** no single scalar is both bounded and size-invariant under class
imbalance. AUPR collapses with size, the AUPR-ratio over-corrects, MCC declines
(though more slowly), and significance transforms (z-score / `-log10(p)`) merely
saturate. **AUROC** is the only prevalence-invariant summary shown, though it is
comparatively insensitive to imbalance (Panel B). The practical guidance:
compare AUPR across studies only at a matched negative-sampling ratio, and
report the honest full-edge-set value rather than a balanced-subsample one
(Panels E–F).

Per-panel data are also written as `grn_reliability_panel*.tsv` for independent
re-plotting.
