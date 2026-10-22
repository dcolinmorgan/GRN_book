#!/usr/bin/env python3
"""Figures for the 50-pair NestBoot SCENIC+ analysis."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SRC = Path("results_summary")
DEST = SRC / "plots"
DEST.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 140,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(DEST / f"{name}.png")
    fig.savefig(DEST / f"{name}.pdf")
    plt.close(fig)


def paired_sizes():
    per = pd.read_csv(SRC / "per_run_counts.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.scatter(per["null_triplets"], per["real_triplets"], s=28, c="#1f4e79", zorder=3)
    lim = max(per["real_triplets"].max(), per["null_triplets"].max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#b0b0b0", lw=1, zorder=1)
    ax.set_xlim(0, per["null_triplets"].max() * 1.15)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Null triplets")
    ax.set_ylabel("Real triplets")
    ax.set_title("Each point is one paired bootstrap")
    save(fig, "triplets_real_vs_null")

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    order = np.argsort(per["fold_triplets"].to_numpy())
    ax.bar(np.arange(len(per)), per["fold_triplets"].to_numpy()[order], color="#1f4e79")
    ax.axhline(per["fold_triplets"].median(), color="#c47b2b", lw=1.2, label=f"median {per['fold_triplets'].median():.1f}×")
    ax.set_xlabel("Runs, ordered by fold")
    ax.set_ylabel("Real / null triplets")
    ax.legend(frameon=False)
    save(fig, "triplet_fold")


def link_support():
    links = pd.read_csv(SRC / "tf_gene_nestboot.tsv", sep="\t", usecols=["real_n", "null_n", "q_mcnemar"])
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    bins = np.arange(0, 52) - 0.5
    ax.hist(links.loc[links.real_n > 0, "real_n"], bins=bins, color="#1f4e79", alpha=0.85, label="Real")
    ax.hist(links.loc[links.null_n > 0, "null_n"], bins=bins, color="#c47b2b", alpha=0.75, label="Null")
    ax.set_xlabel("Runs in which the TF–gene link appears")
    ax.set_ylabel("Links")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    ax.set_title("Recurrence across 50 runs")
    save(fig, "link_recurrence")

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    sig = links["q_mcnemar"] <= 0.05
    ax.scatter(links.loc[~sig, "null_n"], links.loc[~sig, "real_n"], s=6, c="#c8c8c8", rasterized=True, label="q > 0.05")
    ax.scatter(links.loc[sig, "null_n"], links.loc[sig, "real_n"], s=8, c="#1f4e79", rasterized=True, label="q ≤ 0.05")
    ax.set_xlabel("Null runs")
    ax.set_ylabel("Real runs")
    ax.set_xlim(-1, 50)
    ax.set_ylim(-1, 51)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("TF–gene support")
    save(fig, "link_support_scatter")


def top_tfs():
    tf = pd.read_csv(SRC / "tf_significant_targets.tsv", sep="\t").head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    ax.barh(tf["TF"], tf["n_genes"], color="#1f4e79")
    ax.set_xlabel("Target genes with q ≤ 0.05 and real support ≥ 25/50")
    ax.set_title("Stable real-enriched targets")
    save(fig, "top_tfs")


def eregulons():
    eg = pd.read_csv(SRC / "eregulon_nestboot.tsv", sep="\t")
    sig = eg[eg["q_mcnemar"] <= 0.05].sort_values(["real_n", "mean_genes_real"], ascending=[False, False]).head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    y = np.arange(len(sig))
    ax.barh(y - 0.18, sig["real_n"], height=0.36, color="#1f4e79", label="Real")
    ax.barh(y + 0.18, sig["null_n"], height=0.36, color="#c47b2b", label="Null")
    ax.set_yticks(y)
    ax.set_yticklabels(sig["eRegulon"])
    ax.set_xlabel("Runs out of 50")
    ax.set_xlim(0, 52)
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Twenty eRegulons with highest real support, q ≤ 0.05")
    save(fig, "eregulon_support")


if __name__ == "__main__":
    paired_sizes()
    link_support()
    top_tfs()
    eregulons()
    print("wrote", DEST)
