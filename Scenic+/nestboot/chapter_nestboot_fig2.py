#!/usr/bin/env python3
"""Seed support, overlap, and frequency for a SCENIC+ GRNBoost2 benchmark.

Paste into the SCENIC+ chapter. The 50 measured runs differ only in the
GRNBoost2 seed; cells are not resampled. The 50 gene-shuffled runs are the
control. A link's support is the fraction of seeds that contain it, and the
consensus is every measured link whose support exceeds the cutoff set by the
shuffled runs. This is not NestBoot: there is no experiment resampling and no
outer loop.

Expected layout, one directory per run:

    results/real_0/tf_to_gene.tsv
    results/null_0/tf_to_gene.tsv
    ...
    results/real_49/tf_to_gene.tsv
    results/null_49/tf_to_gene.tsv

Each table has columns TF, target, importance. Self-links are dropped because
a TF predicting itself is unchanged by a per-gene shuffle.
"""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS = Path("results")
N_RUNS = 50
# Degree cuts shown in the 2x2. "all" keeps every non-self edge.
PANELS = (
    ("3 TFs per gene", "gene", 3),
    ("5 TFs per gene", "gene", 5),
    ("3 targets per TF", "tf", 3),
    ("5 targets per TF", "tf", 5),
    ("All edges", "all", None),
)


# ---------------------------------------------------------------------------
# 1. Read one GRNBoost2 table and drop TF -> itself.
# ---------------------------------------------------------------------------
def read_edges(path):
    df = pd.read_csv(path, sep="\t", usecols=["TF", "target", "importance"])
    return df.loc[df["TF"].to_numpy() != df["target"].to_numpy()]


# ---------------------------------------------------------------------------
# 2. Sparsity. Keep the strongest n edges of each gene, or of each TF.
#    "all" is the SCENIC+ TF-to-gene table with no degree cut.
#    Within a group, higher importance is the stronger edge.
# ---------------------------------------------------------------------------
def select_edges(df, side, degree):
    if side == "all":
        return df
    group, other = ("target", "TF") if side == "gene" else ("TF", "target")
    parts = []
    for _, block in df.groupby(group, sort=False):
        parts.append(block.nlargest(degree, "importance"))
    return pd.concat(parts, ignore_index=True)


def edge_keys(df):
    return set(zip(df["TF"].to_numpy(), df["target"].to_numpy()))


# ---------------------------------------------------------------------------
# 3. Seed support. For each link, the fraction of the 50 runs that contain it.
#    Returned as a dict {(TF, target): number of runs}.
# ---------------------------------------------------------------------------
def support_counts(side, degree, kind):
    counts = {}
    for i in range(N_RUNS):
        df = read_edges(RESULTS / f"{kind}_{i}" / "tf_to_gene.tsv")
        for key in edge_keys(select_edges(df, side, degree)):
            counts[key] = counts.get(key, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# 4. Cutoff from the shuffled control (Equation 1 of Morgan et al. 2019).
#    FDR(B) is the number of shuffled links with support > B divided by the
#    number of measured links with support > B.
#    The cutoff is the smallest B at which that ratio is at most 0.05.
#    Support > k/50 means the link is in at least k + 1 runs.
# ---------------------------------------------------------------------------
def equation1(measured, shuffled):
    rows = []
    cutoff = None
    for k in range(N_RUNS):
        n_measured = sum(v > k for v in measured.values())
        n_shuffled = sum(v > k for v in shuffled.values())
        fdr = (n_shuffled / n_measured) if n_measured else np.nan
        if cutoff is None and n_measured and fdr <= 0.05:
            cutoff = k / N_RUNS
        rows.append((k / N_RUNS, n_measured, n_shuffled, fdr))
    return rows, cutoff


# ---------------------------------------------------------------------------
# 5. Frequency and overlap, the two y-axes of Figure 2.
#    Frequency at support s is the fraction of links present in exactly
#    s * 50 runs.
#    Overlap at support s is the mean pairwise Jaccard of the 50 run
#    networks after discarding links whose support is below s.
#    Jaccard is |intersection| / |union|.
# ---------------------------------------------------------------------------
def frequency(counts):
    hist = np.bincount(list(counts.values()), minlength=N_RUNS + 1).astype(float)
    hist[0] = 0.0
    total = hist.sum()
    return hist / total if total else hist


def overlap_curve(run_sets, counts):
    """run_sets[i] is the edge set of run i. counts is the pooled support."""
    overlap = np.zeros(N_RUNS + 1)
    for m in range(1, N_RUNS + 1):
        kept = {key for key, n in counts.items() if n >= m}
        sizes = []
        intersections = []
        restricted = [edges & kept for edges in run_sets]
        for a, b in combinations(range(N_RUNS), 2):
            union = len(restricted[a] | restricted[b])
            if union == 0:
                continue
            intersections.append(len(restricted[a] & restricted[b]) / union)
        overlap[m] = float(np.mean(intersections)) if intersections else 0.0
    return overlap


def run_sets(side, degree, kind):
    sets = []
    for i in range(N_RUNS):
        df = read_edges(RESULTS / f"{kind}_{i}" / "tf_to_gene.tsv")
        sets.append(edge_keys(select_edges(df, side, degree)))
    return sets


# ---------------------------------------------------------------------------
# 6. One panel. Frequency only. Each point is the raw fraction of links
#    whose support is exactly that multiple of 1/50. Nothing is smoothed.
#    The line only joins those 50 points. The y-axis runs from 0 to 0.2.
#    The one-run peak sits above that frame. The black line is the
#    cutoff from step 4.
# ---------------------------------------------------------------------------
def draw(ax, title, measured_counts, shuffled_counts, measured_sets, shuffled_sets):
    _, cutoff = equation1(measured_counts, shuffled_counts)
    support = np.arange(1, N_RUNS + 1) / N_RUNS
    measured_frequency = frequency(measured_counts)
    shuffled_frequency = frequency(shuffled_counts)

    ax.plot(support, measured_frequency[1:], color="#2166ac", lw=1.6, marker="o", ms=3, label="Measured")
    ax.plot(support, shuffled_frequency[1:], color="#b2182b", lw=1.6, marker="^", ms=3, label="Shuffled")
    ax.axvline(cutoff, color="black", lw=1.0, label="FDR 0.05")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 0.2)
    ax.set_xlabel("Support")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, ls="--", alpha=0.3)
    return ax


def main():
    # The 2x2 is degree 3 and degree 5, gene side on top and TF side below.
    # The unthresholded table is a separate full-size figure.
    grid = [p for p in PANELS if p[2] in (3, 5)]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.2))
    handles = None
    for ax, (title, side, degree) in zip(axes.ravel(), grid):
        measured = support_counts(side, degree, "real")
        shuffled = support_counts(side, degree, "null")
        measured_sets = run_sets(side, degree, "real")
        shuffled_sets = run_sets(side, degree, "null")
        ax = draw(ax, title, measured, shuffled, measured_sets, shuffled_sets)
        handles = ax.lines
    fig.legend(handles, [h.get_label() for h in handles], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig("fig2_degree_3_and_5.png")

    title, side, degree = PANELS[-1]
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    measured = support_counts(side, degree, "real")
    shuffled = support_counts(side, degree, "null")
    ax = draw(
        ax, title, measured, shuffled, run_sets(side, degree, "real"), run_sets(side, degree, "null")
    )
    handles = ax.lines
    fig.legend(handles, [h.get_label() for h in handles], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig("fig2_all_edges.png")


if __name__ == "__main__":
    main()
