"""
Figure for "How reliable is GRN inference?" chapter.

Illustrates how class imbalance grows with GRN size and why this makes AUPR
(but not AUROC) collapse, and how negative subsampling can inflate a poor AUPR
back to an impressive-looking value — the same scorer, very different numbers.

Core setup: a synthetic scorer whose *discriminative quality is held fixed*
across sizes (true-edge scores ~ N(d,1), non-edge scores ~ N(0,1), same d).
Only the class ratio changes, because in an N-gene network the number of true
edges grows ~N (fixed out-degree k) while possible edges grow ~N^2, so
prevalence ~ k/N falls as the network grows.

Panels:
  A. Prevalence (true-edge fraction) vs. network size — the imbalance driver.
  B. ROC curves for 50/500/5000 genes — nearly identical (deceptively good).
  C. PR curves for the same sizes — collapse as size grows.
  D. How negative subsampling inflates a genuinely poor AUPR: the SAME 5000-gene
     scorer reports AUPR from ~0.02 (full edge set) to ~0.85 (balanced 1:1),
     with the model unchanged.
"""
import csv

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score

rng = np.random.default_rng(42)

# ── Model parameters ──
K_OUTDEGREE = 6      # avg true regulatory edges per gene (biologically ~3-8)
SEPARATION = 1.5     # scorer quality: TP ~ N(SEPARATION,1), FP ~ N(0,1). FIXED across sizes.
SIZES = [50, 500, 5000]
SIZE_COLORS = {50: "#4C72B0", 500: "#DD8452", 5000: "#C44E52"}

# For 5000 genes there are ~25M possible edges; subsample negatives for tractable
# curve computation. This is itself a lesson about "balanced sampling" — but here
# we keep the TRUE prevalence in the metrics by weighting, so the curves are honest.
MAX_NEG_FOR_CURVE = 200_000


def aupr_ratio(prev, aupr):
    """
    AUPR-ratio = observed AUPR / expected-random AUPR. The expected-random AUPR
    equals the prevalence (the area under the flat random-baseline PR line), so
    this is AUPR / prevalence — "how many times better than random."

    NOTE: this is shown for illustration only. It is NOT size-invariant: it
    over-corrects and RISES with network size, because a genuinely
    discriminative scorer's AUPR decays sub-linearly in prevalence. As Erik
    Sonnhammer noted and as DREAM5 (Marbach et al., Nat Methods 2012) implies,
    no AUPR-derived scalar gives a bounded, size-invariant quality number —
    significance transforms (permutation z-score / -log10(p)) instead grow
    astronomically with sample size. AUROC is the only prevalence-invariant
    summary here.
    """
    return aupr / prev


def make_scored_network(N, seed):
    """Return (y_true, scores, prevalence) for an N-gene network with a fixed-quality scorer."""
    r = np.random.default_rng(seed)
    n_possible = N * N - N          # exclude self-loops
    n_pos = K_OUTDEGREE * N
    n_neg = n_possible - n_pos
    prevalence = n_pos / n_possible

    pos_scores = r.normal(SEPARATION, 1.0, n_pos)
    neg_scores = r.normal(0.0, 1.0, n_neg)
    return n_pos, n_neg, pos_scores, neg_scores, prevalence


def metrics_full(n_pos, n_neg, pos_scores, neg_scores):
    """Exact AUROC/AUPR on the full (possibly huge) edge set."""
    y = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
    s = np.concatenate([pos_scores, neg_scores])
    return roc_auc_score(y, s), average_precision_score(y, s)


def mcc_topk(n_pos, n_neg, pos_scores, neg_scores):
    """
    Matthews correlation coefficient (MCC) at the top-k operating point, where
    k = n_pos (predict the k highest-scoring edges as positive — the natural
    threshold for a ranked GRN prediction). MCC is bounded in [-1, 1] and folds
    all four confusion-matrix cells (TP, FP, TN, FN) into one balanced score:

        MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))

    Because it uses TN as well as TP, MCC is far less sensitive to the extreme
    negative inflation that collapses AUPR under class imbalance.
    """
    from math import sqrt
    k = n_pos
    # A score threshold that selects the top-k edges. With TP scores ~N(1.5,1)
    # and FP scores ~N(0,1), rank all scores and call the top k "predicted +".
    all_scores = np.concatenate([pos_scores, neg_scores])
    thresh = np.partition(all_scores, -k)[-k]  # k-th largest score
    tp = int(np.sum(pos_scores >= thresh))
    fp = int(np.sum(neg_scores >= thresh))
    fn = n_pos - tp
    tn = n_neg - fp
    denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def curves_subsampled(n_pos, n_neg, pos_scores, neg_scores, prevalence, seed):
    """
    ROC & PR curves. Subsample negatives if huge, then reweight so the curve
    reflects the TRUE prevalence (precision is corrected back to the real ratio).
    """
    r = np.random.default_rng(seed)
    if n_neg > MAX_NEG_FOR_CURVE:
        neg_scores = r.choice(neg_scores, MAX_NEG_FOR_CURVE, replace=False)
        neg_weight = n_neg / MAX_NEG_FOR_CURVE   # each sampled negative stands for this many
    else:
        neg_weight = 1.0

    y = np.concatenate([np.ones(len(pos_scores)), np.zeros(len(neg_scores))])
    s = np.concatenate([pos_scores, neg_scores])
    w = np.concatenate([np.ones(len(pos_scores)), np.full(len(neg_scores), neg_weight)])

    fpr, tpr, _ = roc_curve(y, s, sample_weight=w)
    prec, rec, _ = precision_recall_curve(y, s, sample_weight=w)
    return (fpr, tpr), (prec, rec)


# ── Compute everything ──
results = {}
for N in SIZES:
    n_pos, n_neg, pos_s, neg_s, prev = make_scored_network(N, seed=100 + N)
    auroc, aupr = metrics_full(n_pos, n_neg, pos_s, neg_s)
    mcc = mcc_topk(n_pos, n_neg, pos_s, neg_s)
    roc_c, pr_c = curves_subsampled(n_pos, n_neg, pos_s, neg_s, prev, seed=200 + N)
    results[N] = dict(prev=prev, auroc=auroc, aupr=aupr, mcc=mcc,
                      aupr_ratio=aupr / prev, roc=roc_c, pr=pr_c)
    print(f"N={N:5d}  prev={prev:.4f}  AUROC={auroc:.3f}  "
          f"AUPR={aupr:.3f}  MCC={mcc:.3f}  AUPR/prev={aupr / prev:.1f}")

# ── Plot ──
fig, axes = plt.subplots(3, 2, figsize=(13, 16))
axA, axB, axC, axD, axE, axF = axes.ravel()

# Panel A: prevalence vs size (the imbalance driver) — sweep a fine size range
size_sweep = np.array([50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000])
prev_sweep = (K_OUTDEGREE * size_sweep) / (size_sweep ** 2 - size_sweep)
axA.loglog(size_sweep, prev_sweep, "-", color="#333333", lw=2, zorder=1)
for N in SIZES:
    axA.plot(N, results[N]["prev"], "o", color=SIZE_COLORS[N], markersize=13,
             zorder=5, label=f"{N} genes")
axA.set_xlabel("Network size (genes)", fontsize=13)
axA.set_ylabel("Prevalence = true edges / possible edges", fontsize=13)
axA.set_title("A. Class imbalance grows with network size", fontsize=13)
axA.legend(fontsize=11, loc="upper right")
axA.grid(True, which="both", alpha=0.3)
axA.annotate("100 genes\n0.03", (100, (K_OUTDEGREE * 100) / (100**2 - 100)),
             xytext=(15, 25), textcoords="offset points", fontsize=10,
             arrowprops=dict(arrowstyle="->", color="gray"))
axA.annotate("1000 genes\n0.003 (10× worse)",
             (1000, (K_OUTDEGREE * 1000) / (1000**2 - 1000)),
             xytext=(15, 25), textcoords="offset points", fontsize=10,
             arrowprops=dict(arrowstyle="->", color="gray"))

# Panel B: ROC curves — nearly identical (deceptively good under imbalance)
for N in SIZES:
    fpr, tpr = results[N]["roc"]
    axB.plot(fpr, tpr, color=SIZE_COLORS[N], lw=2.2,
             label=f"{N} genes (AUROC={results[N]['auroc']:.3f})")
axB.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="random")
axB.set_xlabel("False positive rate", fontsize=13)
axB.set_ylabel("True positive rate", fontsize=13)
axB.set_title("B. ROC curves — nearly identical across sizes", fontsize=13)
axB.legend(fontsize=11, loc="lower right")
axB.grid(True, alpha=0.3)

# Panel C: PR curves — collapse as size grows
for N in SIZES:
    prec, rec = results[N]["pr"]
    axC.plot(rec, prec, color=SIZE_COLORS[N], lw=2.2,
             label=f"{N} genes (AUPR={results[N]['aupr']:.3f})")
    axC.axhline(results[N]["prev"], color=SIZE_COLORS[N], ls="--", lw=1.3, alpha=0.9)
axC.set_xlabel("Recall", fontsize=13)
axC.set_ylabel("Precision", fontsize=13)
axC.set_title("C. PR curves — collapse as size grows", fontsize=13)
axC.legend(fontsize=11, loc="upper right")
axC.set_yscale("log")
axC.grid(True, which="both", alpha=0.3)

# ── Panel D data: how negative subsampling inflates AUPR (same 5000-gene scorer) ──
# Keep all true edges; subsample negatives to a range of neg:pos ratios, from
# balanced (1:1) up to the full edge set. The scorer never changes.
D_N = 5000
_dr = np.random.default_rng(303)
n_pos_D = K_OUTDEGREE * D_N
n_neg_full_D = D_N * D_N - D_N - n_pos_D
pos_D = _dr.normal(SEPARATION, 1.0, n_pos_D)
neg_full_D = _dr.normal(0.0, 1.0, n_neg_full_D)
full_ratio_D = n_neg_full_D / n_pos_D

sampling_ratios = [1, 2, 5, 10, 25, 50, 100, 250, 500, full_ratio_D]
D_rows = []  # (neg_per_pos, prevalence, reported_AUPR, reported_AUROC, aupr_ratio)
for ratio in sampling_ratios:
    n_keep = min(int(round(ratio * n_pos_D)), n_neg_full_D)
    neg_k = _dr.choice(neg_full_D, n_keep, replace=False)
    y = np.concatenate([np.ones(n_pos_D), np.zeros(len(neg_k))])
    s = np.concatenate([pos_D, neg_k])
    aupr_r = average_precision_score(y, s)
    auroc_r = roc_auc_score(y, s)
    prev_r = n_pos_D / (n_pos_D + len(neg_k))
    # AUPR-ratio = actual AUPR / expected-random AUPR; expected-random AUPR = area
    # under the flat baseline at precision=prevalence over recall[0,1] = prevalence
    aupr_ratio_r = aupr_r / prev_r
    D_rows.append((len(neg_k) / n_pos_D, prev_r, aupr_r, auroc_r, aupr_ratio_r))
    print(f"neg:pos={len(neg_k) / n_pos_D:8.1f}:1  prev={prev_r:.4f}  "
          f"AUPR={aupr_r:.3f}  AUROC={auroc_r:.3f}  AUPR/prev={aupr_ratio_r:.1f}")

# Panel D: SAME SCORER, DIFFERENT VERDICTS. A fixed-quality scorer is evaluated on
# networks of increasing size. AUROC (green) stays flat; MCC @ top-k (purple) is more
# balanced but still declines; AUPR (red) collapses with prevalence. The AUPR-ratio
# (yellow dashed, right axis) is the "unfair fix" we critique: dividing AUPR by its
# random baseline OVER-corrects, so it RISES with size instead of staying flat — it is
# not size-invariant either. Only AUROC is genuinely prevalence-invariant here.
d_sizes = np.array(SIZES, dtype=float)
d_auroc_s = np.array([results[N]["auroc"] for N in SIZES])
d_aupr_s = np.array([results[N]["aupr"] for N in SIZES])
d_mcc_s = np.array([results[N]["mcc"] for N in SIZES])
d_ratio_s = np.array([results[N]["aupr_ratio"] for N in SIZES])

l_auroc, = axD.semilogx(d_sizes, d_auroc_s, "-o", color="#55A868", lw=2.2,
                        markersize=9, label="AUROC (prevalence-invariant)")
l_mcc, = axD.semilogx(d_sizes, d_mcc_s, "-s", color="#8172B3", lw=2.2,
                      markersize=9, label="MCC @ top-k (imbalance-robust)")
l_aupr, = axD.semilogx(d_sizes, d_aupr_s, "-o", color="#C44E52", lw=2.2,
                       markersize=9, label="AUPR (collapses with size)")
axD.set_xlabel("Network size (genes)", fontsize=13)
axD.set_ylabel("AUROC / AUPR", fontsize=13)
axD.set_ylim(0, 1)
axD.set_xticks(SIZES)
axD.set_xticklabels([str(N) for N in SIZES])
axD.set_title("D. Same scorer, different verdicts", fontsize=13)
axD.grid(True, which="both", alpha=0.3)

# Secondary axis: AUPR-ratio (AUPR / prevalence) — the "unfair" normalization we
# critique. It over-corrects and RISES with size, so it is not size-invariant.
axD2 = axD.twinx()
l_ratio, = axD2.semilogx(d_sizes, d_ratio_s, "--D", color="#B0A030", lw=2.0,
                         markersize=8, alpha=0.9,
                         label="AUPR-ratio (AUPR / prevalence, over-corrects)")
axD2.set_ylabel("AUPR-ratio (AUPR / prevalence)", fontsize=12, color="#B0A030")
axD2.tick_params(axis="y", labelcolor="#B0A030")
axD2.set_ylim(0, d_ratio_s.max() * 1.35)

axD.legend(handles=[l_auroc, l_mcc, l_aupr, l_ratio], fontsize=9, loc="center left")

# ── Panel E: negative subsampling inflates a poor AUPR — for ALL three sizes ──
# Repeat the Panel-D subsampling sweep for each network size (same fixed-quality
# scorer), showing reported AUPR vs neg:pos ratio per size.
E_rows = {}  # N -> list of (neg_per_pos, prevalence, reported_AUPR)
for N in SIZES:
    _er = np.random.default_rng(404 + N)
    n_pos_E = K_OUTDEGREE * N
    n_neg_full_E = N * N - N - n_pos_E
    pos_E = _er.normal(SEPARATION, 1.0, n_pos_E)
    neg_full_E = _er.normal(0.0, 1.0, n_neg_full_E)
    full_ratio_E = n_neg_full_E / n_pos_E
    ratios_E = [1, 2, 5, 10, 25, 50, 100, 250, 500, full_ratio_E]
    rows = []
    for ratio in ratios_E:
        n_keep = min(int(round(ratio * n_pos_E)), n_neg_full_E)
        neg_k = _er.choice(neg_full_E, n_keep, replace=False)
        y = np.concatenate([np.ones(n_pos_E), np.zeros(len(neg_k))])
        s = np.concatenate([pos_E, neg_k])
        aupr_e = average_precision_score(y, s)
        prev_e = n_pos_E / (n_pos_E + len(neg_k))
        rows.append((len(neg_k) / n_pos_E, prev_e, aupr_e))
    E_rows[N] = rows
    e_ratio = np.array([r[0] for r in rows])
    e_aupr = np.array([r[2] for r in rows])
    e_prev = rows[-1][1]  # prevalence at the full edge set = random-AUPR baseline
    axE.semilogx(e_ratio, e_aupr, "-o", color=SIZE_COLORS[N], lw=2.2,
                 markersize=7, label=f"{N} genes")
    axE.axhline(e_prev, color=SIZE_COLORS[N], ls="--", lw=1.3, alpha=0.9)

axE.set_xlabel("Negatives kept per positive (neg:pos)", fontsize=13)
axE.set_ylabel("Reported AUPR", fontsize=13)
axE.set_yscale("log")
axE.set_title("E. Subsampling inflates a poor AUPR (all sizes)", fontsize=13)
axE.legend(fontsize=11, loc="upper right")
axE.grid(True, which="both", alpha=0.3)

# ── Panel F: same sweep as E, on a LINEAR y-axis ──
# Identical data to Panel E (reported AUPR vs neg:pos ratio, per size), but linear
# instead of log y. On this scale the honest, high-negative end collapses to nearly
# zero while the balanced (1:1) end stays high — showing in true proportion how much
# subsampling inflates AUPR. Prevalence baselines (dashed) sit against the x-axis.
for N in SIZES:
    rows = E_rows[N]
    e_ratio = np.array([r[0] for r in rows])
    e_aupr = np.array([r[2] for r in rows])
    e_prev = rows[-1][1]
    axF.semilogx(e_ratio, e_aupr, "-o", color=SIZE_COLORS[N], lw=2.2,
                 markersize=7, label=f"{N} genes")
    axF.axhline(e_prev, color=SIZE_COLORS[N], ls="--", lw=1.3, alpha=0.9)

axF.set_xlabel("Negatives kept per positive (neg:pos)", fontsize=13)
axF.set_ylabel("Reported AUPR (linear scale)", fontsize=13)
axF.set_ylim(0, 1)
axF.set_title("F. Same as E, linear scale (collapse in true proportion)", fontsize=13)
axF.legend(fontsize=11, loc="upper right")
axF.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("grn_reliability_metrics.png", dpi=150)
plt.savefig("grn_reliability_metrics.pdf")
print("Saved: grn_reliability_metrics.png, grn_reliability_metrics.pdf")

# ── Emit the underlying summary data (one row per evaluated size) ──
with open("grn_reliability_data.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["n_genes", "out_degree_k", "separation_d", "possible_edges",
                "true_edges", "prevalence", "AUROC", "AUPR", "AUPR_ratio"])
    for N in SIZES:
        n_possible = N * N - N
        n_pos = K_OUTDEGREE * N
        r = results[N]
        w.writerow([N, K_OUTDEGREE, SEPARATION, n_possible, n_pos,
                    f"{r['prev']:.6f}", f"{r['auroc']:.4f}",
                    f"{r['aupr']:.4f}", f"{r['aupr_ratio']:.3f}"])
print("Saved: grn_reliability_data.csv")


# ── Per-subfigure data tables (TSV, copy-paste into Google Sheets) ──
def _write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"Saved: {path}")


# Panel A: prevalence vs size (the full swept curve)
_write_tsv(
    "grn_reliability_panelA.tsv",
    ["n_genes", "prevalence"],
    [[int(n), f"{p:.6f}"] for n, p in zip(size_sweep, prev_sweep)],
)

# Panel B: FULL ROC curve points, one file per size (point counts differ per size)
for N in SIZES:
    fpr, tpr = results[N]["roc"]
    _write_tsv(
        f"grn_reliability_panelB_{N}genes.tsv",
        ["FPR", "TPR"],
        [[f"{a:.6f}", f"{b:.6f}"] for a, b in zip(fpr, tpr)],
    )

# Panel C: FULL PR curve points, one file per size
for N in SIZES:
    prec, rec = results[N]["pr"]
    _write_tsv(
        f"grn_reliability_panelC_{N}genes.tsv",
        ["recall", "precision"],
        [[f"{r_:.6f}", f"{p_:.6f}"] for r_, p_ in zip(rec, prec)],
    )

# Panel C baselines (dashed lines = prevalence per size)
_write_tsv(
    "grn_reliability_panelC_baselines.tsv",
    ["n_genes", "random_baseline_precision"],
    [[N, f"{results[N]['prev']:.6f}"] for N in SIZES],
)

# Panel D: same-scorer verdicts across network sizes (AUROC, AUPR, AUPR-ratio, AUPR z-score)
_write_tsv(
    "grn_reliability_panelD.tsv",
    ["n_genes", "prevalence", "AUROC", "AUPR", "MCC_topk", "AUPR_ratio"],
    [[N, f"{results[N]['prev']:.6f}", f"{results[N]['auroc']:.4f}",
      f"{results[N]['aupr']:.4f}", f"{results[N]['mcc']:.4f}",
      f"{results[N]['aupr_ratio']:.3f}"]
     for N in SIZES],
)

# Panel D (neg:pos subsampling sweep, retained as separate data for reference)
_write_tsv(
    "grn_reliability_panelD_subsampling.tsv",
    ["neg_per_pos", "prevalence", "reported_AUPR", "reported_AUROC", "AUPR_ratio"],
    [[f"{r[0]:.2f}", f"{r[1]:.6f}", f"{r[2]:.4f}", f"{r[3]:.4f}", f"{r[4]:.3f}"] for r in D_rows],
)

# Panel E: reported AUPR vs negative-sampling ratio, one file per size
for N in SIZES:
    _write_tsv(
        f"grn_reliability_panelE_{N}genes.tsv",
        ["neg_per_pos", "prevalence", "reported_AUPR"],
        [[f"{r[0]:.2f}", f"{r[1]:.6f}", f"{r[2]:.4f}"] for r in E_rows[N]],
    )

# Panel F uses the same data as Panel E (reported AUPR vs neg:pos ratio per size),
# only on a linear y-axis; see grn_reliability_panelE_*.tsv.
