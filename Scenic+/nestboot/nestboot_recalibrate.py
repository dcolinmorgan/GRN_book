#!/usr/bin/env python3
"""NestBoot Eq. 1 at 3 links per node, using the 50 seeds as one inner loop.

Self-links are removed before the degree cut. A TF predicting itself is a
perfect copy of the label and is unchanged by the gene shuffle. The 2019
paper does not treat that edge as a regulatory link.

Q = 50 GRNBoost2 seeds. There is no outer loop, so the consensus is every
link with support > B. Support > k/50 means the link is in at least k+1 seeds.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
PLOTS = DEST / "plots"
N = 50
DEGREE = 3


def encode(values, vocab):
    out = np.empty(len(values), dtype=np.int32)
    for i, v in enumerate(values):
        code = vocab.get(v)
        if code is None:
            code = len(vocab)
            vocab[v] = code
        out[i] = code
    return out


def top_per_group(group, other, imp, n):
    order = np.lexsort((other, -imp, group))
    g = group[order]
    o = other[order]
    change = np.flatnonzero(g[1:] != g[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(g)]
    take = np.concatenate([np.arange(s, min(s + n, e)) for s, e in zip(starts, ends)])
    return g[take], o[take]


def pack(tf, gene):
    return tf.astype(np.int64) * 100_000 + gene.astype(np.int64)


def hist(arrays):
    stacked = np.concatenate([np.unique(a) for a in arrays])
    stacked.sort()
    change = np.flatnonzero(stacked[1:] != stacked[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, stacked.size]
    support = ends - starts
    return np.bincount(support, minlength=N + 1).astype(np.int64)


def eq1(real_h, null_h):
    rows = []
    hit = None
    for k in range(0, N):
        n_real = int(real_h[k + 1 :].sum())
        n_null = int(null_h[k + 1 :].sum())
        fdr = (n_null / n_real) if n_real else np.nan
        row = {
            "B_runs": k,
            "B": k / N,
            "support_gt": k + 1,
            "measured": n_real,
            "shuffled": n_null,
            "fdr": fdr,
        }
        rows.append(row)
        if hit is None and n_real and fdr <= 0.05:
            hit = row
    return rows, hit


def main():
    DEST.mkdir(exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    vocab_tf, vocab_gene = {}, {}
    gene_real, gene_null = [], []
    tf_real, tf_null = [], []
    n_gene = {"real": [], "null": []}
    n_tf = {"real": [], "null": []}
    for i in range(N):
        for kind, gstore, tstore in (
            ("real", gene_real, tf_real),
            ("null", gene_null, tf_null),
        ):
            df = pd.read_csv(
                OUT / f"{kind}_{i}" / "tf_to_gene.tsv",
                sep="\t",
                usecols=["TF", "target", "importance"],
            )
            df = df[df["TF"].to_numpy() != df["target"].to_numpy()]
            tf = encode(df["TF"].to_numpy(), vocab_tf)
            gene = encode(df["target"].to_numpy(), vocab_gene)
            imp = df["importance"].to_numpy()
            gg, tt = top_per_group(gene, tf, imp, DEGREE)
            gstore.append(pack(tt, gg))
            n_gene[kind].append(len(gg))
            tt2, gg2 = top_per_group(tf, gene, imp, DEGREE)
            tstore.append(pack(tt2, gg2))
            n_tf[kind].append(len(gg2))
        print(f"run {i}", flush=True)

    blocks = []
    for name, real_a, null_a, counts in (
        ("3 TFs per gene, self-links removed", gene_real, gene_null, n_gene),
        ("3 targets per TF, self-links removed", tf_real, tf_null, n_tf),
    ):
        rh, nh = hist(real_a), hist(null_a)
        rows, hit = eq1(rh, nh)
        blocks.append((name, counts, rh, nh, rows, hit))
        print(name, "edges", int(np.median(counts["real"])), "cutoff", hit, flush=True)
        pd.DataFrame(rows).to_csv(
            DEST / ("eq1_3per_gene.tsv" if "gene" in name else "eq1_3per_tf.tsv"),
            sep="\t",
            index=False,
        )

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 200,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))
    colors = ["#1f4e79", "#c47b2b"]
    for ax_i, (name, counts, rh, nh, rows, hit) in enumerate(blocks):
        support = np.arange(1, N + 1) / N
        ax = axes[0]
        ax.plot(support, rh[1:], color=colors[ax_i], lw=1.6, label=f"Measured, {name.split(',')[0]}")
        ax.plot(support, nh[1:], color=colors[ax_i], lw=1.2, ls="--", label=f"Shuffled, {name.split(',')[0]}")
        if hit is not None:
            ax.axvline(hit["B"], color=colors[ax_i], ls=":", lw=1)
        ax2 = axes[1]
        ax2.plot([r["B"] for r in rows], [r["fdr"] for r in rows], color=colors[ax_i], lw=1.6, label=name.split(",")[0])
        if hit is not None:
            ax2.scatter([hit["B"]], [hit["fdr"]], color=colors[ax_i], zorder=3)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Bootstrap support")
    axes[0].set_ylabel("Links at this support")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(True, ls="--", alpha=0.3)
    axes[0].set_title("Link frequency")
    axes[1].axhline(0.05, color="gray", ls=":", lw=1)
    axes[1].set_xlabel("Support cutoff B")
    axes[1].set_ylabel("Shuffled links > B  /  measured links > B")
    axes[1].set_ylim(0, 1.35)
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(True, ls="--", alpha=0.3)
    axes[1].set_title("Equation 1")
    fig.tight_layout()
    fig.savefig(PLOTS / "nestboot_eq1_degree3.png")
    fig.savefig(PLOTS / "nestboot_eq1_degree3.pdf")
    print("wrote figure", flush=True)


if __name__ == "__main__":
    main()
