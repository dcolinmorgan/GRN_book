#!/usr/bin/env python3
"""Figure 2 of Morgan et al. 2019, on the 50 GRNBoost2 seeds.

X: bootstrap support.
Left Y: fraction of links with that exact support (measured and shuffled).
Right Y: mean pairwise Jaccard of the run-networks restricted to links whose
support exceeds the cutoff. Jaccard is |intersection| / |union|.
The vertical line is the Equation 1 cutoff, FDR 0.05.
"""

from itertools import combinations
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
    take = np.concatenate([np.arange(s, min(s + n, int(e))) for s, e in zip(starts, ends)])
    return g[take], o[take]


def pack(tf, gene):
    return tf.astype(np.int64) * 100_000 + gene.astype(np.int64)


def add_run(store, keys, run):
    bit = np.uint64(1) << np.uint64(run)
    for key in keys:
        store[int(key)] = store.get(int(key), np.uint64(0)) | bit


def support_and_overlap(store):
    """Exact-support histogram and mean pairwise Jaccard above each cutoff."""
    bits = np.fromiter(store.values(), dtype=np.uint64, count=len(store))
    # popcount
    supp = np.array([int(b.bit_count()) for b in bits], dtype=np.int16)
    hist = np.bincount(supp, minlength=N + 1).astype(np.int64)
    # group indices by support
    order = np.argsort(supp)
    bits = bits[order]
    supp = supp[order]
    # running pair intersections and per-run sizes, from high support down
    pair = np.zeros((N, N), dtype=np.int64)
    size = np.zeros(N, dtype=np.int64)
    overlap = np.full(N + 1, np.nan)  # overlap[m] = Jaccard of links with support >= m
    # walk unique support values descending
    # precompute run index lists
    run_idx = []
    for b in bits:
        run_idx.append(tuple(i for i in range(N) if (int(b) >> i) & 1))
    # process from the end (highest support)
    start = len(bits)
    for m in range(N, 0, -1):
        # add links whose support == m
        lo = int(np.searchsorted(supp, m, side="left"))
        for t in range(lo, start):
            runs = run_idx[t]
            for i in runs:
                size[i] += 1
            for a, b in combinations(runs, 2):
                pair[a, b] += 1
        start = lo
        jacs = []
        for a, b in combinations(range(N), 2):
            inter = pair[a, b]
            union = size[a] + size[b] - inter
            if union == 0:
                continue
            jacs.append(inter / union)
        overlap[m] = float(np.mean(jacs)) if jacs else 0.0
    return hist, overlap


def first_fdr(real_h, null_h):
    for k in range(0, N):
        n_real = int(real_h[k + 1 :].sum())
        n_null = int(null_h[k + 1 :].sum())
        if n_real and n_null / n_real <= 0.05:
            return k / N, n_real, n_null, n_null / n_real
    return None, None, None, None


def load_graphs():
    vocab_tf, vocab_gene = {}, {}
    gene = {"real": {}, "null": {}}
    tf = {"real": {}, "null": {}}
    for i in range(N):
        for kind in ("real", "null"):
            df = pd.read_csv(
                OUT / f"{kind}_{i}" / "tf_to_gene.tsv",
                sep="\t",
                usecols=["TF", "target", "importance"],
            )
            df = df[df["TF"].to_numpy() != df["target"].to_numpy()]
            tfc = encode(df["TF"].to_numpy(), vocab_tf)
            genec = encode(df["target"].to_numpy(), vocab_gene)
            imp = df["importance"].to_numpy()
            gg, tt = top_per_group(genec, tfc, imp, DEGREE)
            add_run(gene[kind], pack(tt, gg), i)
            tt2, gg2 = top_per_group(tfc, genec, imp, DEGREE)
            add_run(tf[kind], pack(tt2, gg2), i)
        print(f"run {i}", flush=True)
    return gene, tf


def draw(ax, title, real_h, null_h, real_ov, null_ov, b_star):
    support = np.arange(1, N + 1) / N
    freq_r = real_h[1:] / real_h[1:].sum()
    freq_n = null_h[1:] / null_h[1:].sum()
    ax.plot(support, real_ov[1:], color="#2166ac", lw=1.4, ls="--", label="Measured overlap")
    ax.plot(support, null_ov[1:], color="#b2182b", lw=1.4, ls="--", label="Shuffled overlap")
    ax.set_xlabel("Bootstrap support")
    ax.set_ylabel("Overlap  |intersection| / |union|")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.0)
    right = ax.twinx()
    right.plot(support, freq_r, color="#2166ac", lw=1.8, label="Measured frequency")
    right.plot(support, freq_n, color="#b2182b", lw=1.8, label="Shuffled frequency")
    right.set_ylabel("Link frequency")
    right.set_ylim(0, 0.02)
    if b_star is not None:
        ax.axvline(b_star, color="black", lw=1.2, label="FDR 0.05 cutoff")
    ax.set_title(title)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = right.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper center")
    ax.grid(True, ls="--", alpha=0.3)


def main():
    PLOTS.mkdir(exist_ok=True)
    gene, tf = load_graphs()
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "figure.dpi": 140,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    rows = []
    for ax, name, store in (
        (axes[0], "3 TFs per gene", gene),
        (axes[1], "3 targets per TF", tf),
    ):
        print(f"overlap {name} measured", flush=True)
        rh, rov = support_and_overlap(store["real"])
        print(f"overlap {name} shuffled", flush=True)
        nh, nov = support_and_overlap(store["null"])
        b_star, n_real, n_null, fdr = first_fdr(rh, nh)
        draw(ax, name, rh, nh, rov, nov, b_star)
        print(name, "B", b_star, "kept", n_real, n_null, "fdr", fdr, flush=True)
        for m in range(1, N + 1):
            rows.append(
                {
                    "graph": name,
                    "support": m / N,
                    "measured_frequency": rh[m] / rh[1:].sum(),
                    "shuffled_frequency": nh[m] / nh[1:].sum(),
                    "measured_overlap": rov[m],
                    "shuffled_overlap": nov[m],
                }
            )
    fig.tight_layout()
    fig.savefig(PLOTS / "fig2_support_frequency_overlap.png")
    fig.savefig(PLOTS / "fig2_support_frequency_overlap.pdf")
    pd.DataFrame(rows).to_csv(DEST / "fig2_support_frequency_overlap.tsv", sep="\t", index=False)
    print("wrote", PLOTS / "fig2_support_frequency_overlap.png", flush=True)


if __name__ == "__main__":
    main()
