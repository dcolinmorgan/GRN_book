#!/usr/bin/env python3
"""Figure 2 panels for 3–5 links per node, split by gene and by TF, plus all edges.

Self-links are removed. A link's support is the fraction of the 50 seeds that
keep it. Overlap at a support x is the mean pairwise Jaccard of the 50 run
networks after dropping links with support below x. Frequency is the fraction
of links whose support is exactly x.
"""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numba
import numpy as np
import pandas as pd

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
PLOTS = DEST / "plots"
N = 50
DEGREES = (3, 4, 5)


def encode(values, vocab):
    out = np.empty(len(values), dtype=np.int32)
    for i, value in enumerate(values):
        code = vocab.get(value)
        if code is None:
            code = len(vocab)
            vocab[value] = code
        out[i] = code
    return out


def top_per_group(group, other, imp, n):
    order = np.lexsort((other, -imp, group))
    g = group[order]
    o = other[order]
    change = np.flatnonzero(g[1:] != g[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(g)]
    take = np.concatenate([np.arange(s, min(s + n, int(e))) for s, e in zip(starts, ends)])
    return o[take], g[take]


def pack(tf, gene):
    return tf.astype(np.int64) * 100_000 + gene.astype(np.int64)


def add_run(store, keys, run):
    bit = np.uint64(1) << np.uint64(run)
    for key in keys:
        store[int(key)] = store.get(int(key), np.uint64(0)) | bit


@numba.njit
def _overlap_from_bits(bits, supp):
    """bits, supp sorted by supp ascending. overlap[m] is Jaccard at support >= m."""
    pair = np.zeros((N, N), dtype=np.int64)
    size = np.zeros(N, dtype=np.int64)
    overlap = np.zeros(N + 1)
    idx = len(bits) - 1
    for m in range(N, 0, -1):
        while idx >= 0 and supp[idx] == m:
            b = bits[idx]
            runs = np.empty(N, dtype=np.int16)
            c = 0
            for i in range(N):
                if (b >> i) & np.uint64(1):
                    runs[c] = i
                    c += 1
            for a in range(c):
                size[runs[a]] += 1
            for a in range(c):
                for b2 in range(a + 1, c):
                    pair[runs[a], runs[b2]] += 1
            idx -= 1
        num = 0.0
        count = 0
        for a in range(N):
            for b2 in range(a + 1, N):
                union = size[a] + size[b2] - pair[a, b2]
                if union == 0:
                    continue
                num += pair[a, b2] / union
                count += 1
        overlap[m] = num / count if count else 0.0
    return overlap


def curves(store):
    bits = np.fromiter(store.values(), dtype=np.uint64, count=len(store))
    supp = np.array([int(b.bit_count()) for b in bits], dtype=np.int16)
    hist = np.bincount(supp, minlength=N + 1).astype(np.int64)
    order = np.argsort(supp, kind="mergesort")
    overlap = _overlap_from_bits(bits[order], supp[order])
    freq = np.zeros(N + 1)
    total = int(hist[1:].sum())
    if total:
        freq[1:] = hist[1:] / total
    return hist, freq, overlap


def first_fdr(hist_real, hist_null):
    for k in range(N):
        n_real = int(hist_real[k + 1 :].sum())
        n_null = int(hist_null[k + 1 :].sum())
        if n_real and n_null / n_real <= 0.05:
            return k / N, n_real, n_null, n_null / n_real
    return None, None, None, None


def draw(ax, title, real, null, b_star):
    support = np.arange(1, N + 1) / N
    ax.plot(support, real[2][1:], color="#2166ac", ls="--", lw=1.3)
    ax.plot(support, null[2][1:], color="#b2182b", ls="--", lw=1.3)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=10)
    ax.grid(True, ls="--", alpha=0.3)
    right = ax.twinx()
    right.plot(support, real[1][1:], color="#2166ac", lw=1.4)
    right.plot(support, null[1][1:], color="#b2182b", lw=1.4)
    right.set_ylim(0, 0.02)
    if b_star is not None:
        ax.axvline(b_star, color="black", lw=1.0)
    return right


def main():
    PLOTS.mkdir(exist_ok=True)
    vocab_tf, vocab_gene = {}, {}
    graphs = {}
    for degree in DEGREES:
        graphs[(f"{degree} TFs per gene", degree)] = {"real": {}, "null": {}}
        graphs[(f"{degree} targets per TF", degree)] = {"real": {}, "null": {}}
    full = {"real": {}, "null": {}}
    # Stable handles: rebuild by degree inside the loop using fresh stores.
    gene = {d: {"real": {}, "null": {}} for d in DEGREES}
    tf = {d: {"real": {}, "null": {}} for d in DEGREES}

    for i in range(N):
        for kind in ("real", "null"):
            df = pd.read_csv(
                OUT / f"{kind}_{i}" / "tf_to_gene.tsv",
                sep="\t",
                usecols=["TF", "target", "importance"],
            )
            keep = df["TF"].to_numpy() != df["target"].to_numpy()
            df = df.loc[keep]
            tfc = encode(df["TF"].to_numpy(), vocab_tf)
            genec = encode(df["target"].to_numpy(), vocab_gene)
            imp = df["importance"].to_numpy()
            add_run(full[kind], pack(tfc, genec), i)
            for d in DEGREES:
                gg, tt = top_per_group(genec, tfc, imp, d)
                add_run(gene[d][kind], pack(tt, gg), i)
                tt2, gg2 = top_per_group(tfc, genec, imp, d)
                add_run(tf[d][kind], pack(tt2, gg2), i)
        print(f"run {i}", flush=True)

    panels = []
    rows = []
    for d in DEGREES:
        for label, bucket in ((f"{d} TFs per gene", gene[d]), (f"{d} targets per TF", tf[d])):
            print("curves", label, flush=True)
            real = curves(bucket["real"])
            null = curves(bucket["null"])
            hit = first_fdr(real[0], null[0])
            panels.append((label, real, null, hit[0]))
            print(label, "B", hit, flush=True)
            for m in range(1, N + 1):
                rows.append(
                    {
                        "graph": label,
                        "support": m / N,
                        "measured_frequency": real[1][m],
                        "shuffled_frequency": null[1][m],
                        "measured_overlap": real[2][m],
                        "shuffled_overlap": null[2][m],
                    }
                )
    print("curves all edges", flush=True)
    real = curves(full["real"])
    null = curves(full["null"])
    hit = first_fdr(real[0], null[0])
    panels.append(("All edges", real, null, hit[0]))
    print("All edges", "B", hit, flush=True)
    for m in range(1, N + 1):
        rows.append(
            {
                "graph": "All edges",
                "support": m / N,
                "measured_frequency": real[1][m],
                "shuffled_frequency": null[1][m],
                "measured_overlap": real[2][m],
                "shuffled_overlap": null[2][m],
            }
        )

    pd.DataFrame(rows).to_csv(DEST / "fig2_degree_and_full.tsv", sep="\t", index=False)

    order = [
        "3 TFs per gene",
        "4 TFs per gene",
        "5 TFs per gene",
        "3 targets per TF",
        "4 targets per TF",
        "5 targets per TF",
        "All edges",
    ]
    by_name = {name: (real, null, b) for name, real, null, b in panels}
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "savefig.dpi": 200, "savefig.bbox": "tight"}
    )
    fig, axes = plt.subplots(2, 4, figsize=(14.5, 7.2))
    axes = axes.ravel()
    rights = []
    for ax, name in zip(axes, order):
        real, null, b = by_name[name]
        rights.append(draw(ax, name, real, null, b))
    axes[-1].axis("off")
    for ax in axes[:4]:
        ax.set_xlabel("Support")
    for ax in axes[4:7]:
        ax.set_xlabel("Support")
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Overlap")
    for right in (rights[3], rights[6]):
        right.set_ylabel("Frequency")
    fig.tight_layout()
    fig.savefig(PLOTS / "fig2_degree_and_full.png")
    fig.savefig(PLOTS / "fig2_degree_and_full.pdf")
    print("wrote", PLOTS / "fig2_degree_and_full.png", flush=True)


if __name__ == "__main__":
    main()
