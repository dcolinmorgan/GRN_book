#!/usr/bin/env python3
"""Nested bootstrap on the 50 existing eRegulon runs.

The 50 seeds are split into R outer groups of Q inner runs (5 x 10 and 10 x 5).
Groups are interleaved, so outer r contains runs r, r+R, r+2R, ...
Within an outer group, a link's support is the fraction of inner runs that contain it.
The link passes that outer group at cutoff B when support > B.
It is kept only when it passes every outer group.
FDR(B) = shuffled links kept / measured links kept.

Five outer groups of 10 are the primary split: support is estimated to 0.1.
Ten groups of 5 are the check: the same 50 runs, a coarser support grid, and a
stricter every-outer rule.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
PLOTS = DEST / "plots"
N = 50
SCHEMES = ((5, 10), (10, 5))


def load_bits():
    real, null = {}, {}
    for i in range(N):
        for kind, store in (("real", real), ("null", null)):
            df = pd.read_csv(
                OUT / f"{kind}_{i}" / "eRegulons.tsv",
                sep="\t",
                usecols=["TF", "Gene"],
            ).drop_duplicates()
            bit = np.uint64(1) << np.uint64(i)
            for tf, gene in zip(df["TF"].to_numpy(), df["Gene"].to_numpy()):
                key = (tf, gene)
                store[key] = store.get(key, np.uint64(0)) | bit
        print(f"loaded run {i}  measured {len(real)} shuffled {len(null)}", flush=True)
    return real, null


def group_masks(n_outer, n_inner):
    masks = []
    for r in range(n_outer):
        runs = [r + n_outer * j for j in range(n_inner)]
        mask = np.uint64(0)
        for i in runs:
            mask |= np.uint64(1) << np.uint64(i)
        masks.append(mask)
    return masks


def passes_all(bits, masks, min_count):
    """True when every outer group has at least min_count inner hits."""
    for mask in masks:
        if int((bits & mask).bit_count()) < min_count:
            return False
    return True


def evaluate(store, masks, n_inner):
    """For each inner-count threshold m = 1..n_inner, links passing every outer."""
    kept = []
    items = list(store.items())
    for m in range(1, n_inner + 1):
        keys = [key for key, bits in items if passes_all(bits, masks, m)]
        kept.append(keys)
    return kept


def main():
    DEST.mkdir(exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    real, null = load_bits()
    rows = []
    curves = {}
    chosen = {}
    for n_outer, n_inner in SCHEMES:
        masks = group_masks(n_outer, n_inner)
        real_kept = evaluate(real, masks, n_inner)
        null_kept = evaluate(null, masks, n_inner)
        label = f"{n_outer}x{n_inner}"
        curves[label] = []
        hit = None
        for m, (rk, nk) in enumerate(zip(real_kept, null_kept), start=1):
            # support > (m-1)/Q  <=>  at least m inner hits
            b = (m - 1) / n_inner
            fdr = (len(nk) / len(rk)) if rk else np.nan
            row = {
                "scheme": label,
                "n_outer": n_outer,
                "n_inner": n_inner,
                "min_inner_hits": m,
                "B": b,
                "measured_kept": len(rk),
                "shuffled_kept": len(nk),
                "fdr": fdr,
            }
            rows.append(row)
            curves[label].append(row)
            if hit is None and rk and fdr <= 0.05:
                hit = row
                chosen[label] = (rk, nk, row)
            print(
                f"{label}  inner>={m}/{n_inner}  B={b:.2f}  "
                f"measured {len(rk)}  shuffled {len(nk)}  fdr {fdr}",
                flush=True,
            )
        if hit is None:
            print(f"{label} never reaches FDR 0.05", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(DEST / "inner_bootstrap_fdr.tsv", sep="\t", index=False)

    # Overlap of the two selected measured sets.
    if "5x10" in chosen and "10x5" in chosen:
        a = set(chosen["5x10"][0])
        b = set(chosen["10x5"][0])
        inter = len(a & b)
        print(
            f"selected-set overlap  {inter}  "
            f"only 5x10 {len(a - b)}  only 10x5 {len(b - a)}  "
            f"jaccard {inter / len(a | b):.4f}",
            flush=True,
        )
        both = sorted(a & b)
        pd.DataFrame(both, columns=["TF", "Gene"]).to_csv(
            DEST / "inner_bootstrap_links_both.tsv", sep="\t", index=False
        )
        pd.DataFrame(sorted(a), columns=["TF", "Gene"]).to_csv(
            DEST / "inner_bootstrap_links_5x10.tsv", sep="\t", index=False
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
    colors = {"5x10": "#1f4e79", "10x5": "#c47b2b"}
    for label, series in curves.items():
        axes[0].plot(
            [r["B"] for r in series],
            [r["fdr"] for r in series],
            marker="o",
            ms=4,
            color=colors[label],
            label=f"{label.replace('x', ' outer × ')} inner",
        )
        axes[1].plot(
            [r["B"] for r in series],
            [r["measured_kept"] for r in series],
            marker="o",
            ms=4,
            color=colors[label],
            label=f"{label} measured",
        )
        axes[1].plot(
            [r["B"] for r in series],
            [r["shuffled_kept"] for r in series],
            marker="^",
            ms=4,
            ls="--",
            color=colors[label],
            label=f"{label} shuffled",
        )
    axes[0].axhline(0.05, color="gray", ls=":", lw=1)
    axes[0].set_xlabel("Inner-support cutoff B")
    axes[0].set_ylabel("Shuffled kept / measured kept")
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(True, ls="--", alpha=0.3)
    axes[0].set_title("FDR after the link passes every outer group")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Inner-support cutoff B")
    axes[1].set_ylabel("Links kept")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(True, ls="--", alpha=0.3)
    axes[1].set_title("Links kept")
    fig.tight_layout()
    fig.savefig(PLOTS / "inner_bootstrap_fdr.png")
    fig.savefig(PLOTS / "inner_bootstrap_fdr.pdf")
    print("wrote", DEST / "inner_bootstrap_fdr.tsv", flush=True)


if __name__ == "__main__":
    main()
