#!/usr/bin/env python3
"""One Figure 2 curve from the union of the two degree-3 graphs.

A link is kept in a run if it is among that gene's 3 strongest TFs or among
that TF's 3 strongest targets. Self-links are dropped first.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from fig2_overlap import (
    DEGREE,
    DEST,
    N,
    OUT,
    PLOTS,
    add_run,
    encode,
    first_fdr,
    pack,
    support_and_overlap,
    top_per_group,
)


def main():
    vocab_tf, vocab_gene = {}, {}
    store = {"real": {}, "null": {}}
    sizes = {"real": [], "null": []}
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
            tt2, gg2 = top_per_group(tfc, genec, imp, DEGREE)
            keys = set(pack(tt, gg).tolist()) | set(pack(tt2, gg2).tolist())
            add_run(store[kind], keys, i)
            sizes[kind].append(len(keys))
        print(f"run {i}  edges {sizes['real'][-1]}", flush=True)

    print("measured links", len(store["real"]), "shuffled", len(store["null"]), flush=True)
    rh, rov = support_and_overlap(store["real"])
    nh, nov = support_and_overlap(store["null"])
    b_star, n_real, n_null, fdr = first_fdr(rh, nh)
    print("cutoff", b_star, n_real, n_null, fdr, flush=True)
    print("edges per run", min(sizes["real"]), int(pd.Series(sizes["real"]).median()), max(sizes["real"]), flush=True)

    support = [m / N for m in range(1, N + 1)]
    rows = pd.DataFrame(
        {
            "support": support,
            "measured_frequency": rh[1:] / rh[1:].sum(),
            "shuffled_frequency": nh[1:] / nh[1:].sum(),
            "measured_overlap": rov[1:],
            "shuffled_overlap": nov[1:],
        }
    )
    rows.to_csv(DEST / "fig2_union.tsv", sep="\t", index=False)

    plt.rcParams.update(
        {"font.size": 11, "axes.spines.top": False, "savefig.dpi": 200, "savefig.bbox": "tight"}
    )
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(rows.support, rows.measured_overlap, color="#2166ac", ls="--", lw=1.6, label="Measured")
    ax.plot(rows.support, rows.shuffled_overlap, color="#b2182b", ls="--", lw=1.6, label="Shuffled")
    ax.set_xlabel("Support")
    ax.set_ylabel("Overlap")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1)
    right = ax.twinx()
    right.plot(rows.support, rows.measured_frequency, color="#2166ac", lw=1.8)
    right.plot(rows.support, rows.shuffled_frequency, color="#b2182b", lw=1.8)
    right.set_ylabel("Frequency")
    right.set_ylim(0, 0.02)
    if b_star is not None:
        ax.axvline(b_star, color="black", lw=1.1)
    ax.legend(frameon=False, loc="center right")
    ax.grid(True, ls="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "fig2_support_frequency_overlap.png")
    fig.savefig(PLOTS / "fig2_support_frequency_overlap.pdf")
    inside = rows[(rows.shuffled_frequency > 0) & (rows.shuffled_frequency <= 0.02)]
    zero = rows[rows.shuffled_frequency == 0]
    print("shuffled freq in (0, 0.02]", inside[["support", "shuffled_frequency"]].to_string(index=False), flush=True)
    print("first zero shuffled freq", float(zero.support.min()) if len(zero) else None, flush=True)


if __name__ == "__main__":
    main()
