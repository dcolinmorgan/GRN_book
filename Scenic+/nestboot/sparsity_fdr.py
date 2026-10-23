#!/usr/bin/env python3
"""NestBoot FDR at 30 GRNBoost2 sparsities.

Each run keeps its k strongest TF–gene edges by GRNBoost2 importance.
k is log-spaced from 2,000 to 2,000,000 (30 values). Density is k divided
by the complete bipartite graph: 824 TFs x 36,601 genes.

FDR(B) = #{shuffled links with support > B} / #{measured links with support > B}.
Support is the number of runs in which the link is among that run's top k.
"""

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
N_RUNS = 50
N_TF = 824
N_GENE = 36601
N_COMPLETE = N_TF * N_GENE
K_MAX = 2_000_000
KS = np.unique(np.round(np.geomspace(2_000, K_MAX, 30)).astype(np.int64))


def encode(values, vocab):
    out = np.empty(len(values), dtype=np.int32)
    for i, v in enumerate(values):
        code = vocab.get(v)
        if code is None:
            code = len(vocab)
            vocab[v] = code
        out[i] = code
    return out


def top_keys(path, vocab_tf, vocab_gene):
    df = pd.read_csv(path, sep="\t", usecols=["TF", "target", "importance"])
    imp = df["importance"].to_numpy()
    if len(imp) < K_MAX:
        raise SystemExit(f"{path} has {len(imp)} edges, fewer than {K_MAX}")
    idx = np.argpartition(imp, -K_MAX)[-K_MAX:]
    idx = idx[np.argsort(imp[idx], kind="mergesort")[::-1]]
    tf = encode(df["TF"].to_numpy()[idx], vocab_tf)
    gene = encode(df["target"].to_numpy()[idx], vocab_gene)
    return tf.astype(np.int64) * 100_000 + gene.astype(np.int64)


def support_hist(keys, k):
    """keys: (n_runs, K_MAX) int64, ranked strongest-first."""
    stacked = np.concatenate([keys[i, :k] for i in range(keys.shape[0])])
    stacked.sort()
    # run-length encode; a link cannot repeat inside one run's top k
    change = np.flatnonzero(stacked[1:] != stacked[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(stacked)]
    support = ends - starts
    hist = np.bincount(support, minlength=N_RUNS + 1)
    return hist


def first_fdr(real_hist, null_hist):
    """Smallest m such that links in >= m runs have FDR <= 0.05.

    Equivalent to Eq. 1 at B = m - 1 (support > B).
    """
    real_tail = np.cumsum(real_hist[::-1])[::-1]
    null_tail = np.cumsum(null_hist[::-1])[::-1]
    for m in range(1, N_RUNS + 1):
        rr, nn = int(real_tail[m]), int(null_tail[m])
        if rr == 0:
            continue
        fdr = nn / rr
        if fdr <= 0.05:
            return m, fdr, rr, nn
    return None, None, None, None


def main():
    if len(KS) != 30:
        raise SystemExit(f"expected 30 sparsities, got {len(KS)}: {KS}")
    DEST.mkdir(exist_ok=True)
    cache = DEST / "sparsity_top_keys.npz"
    if cache.is_file():
        z = np.load(cache)
        real, null = z["real"], z["null"]
        print(f"loaded {cache}", flush=True)
    else:
        vocab_tf, vocab_gene = {}, {}
        real = np.empty((N_RUNS, K_MAX), dtype=np.int64)
        null = np.empty((N_RUNS, K_MAX), dtype=np.int64)
        for i in range(N_RUNS):
            real[i] = top_keys(OUT / f"real_{i}" / "tf_to_gene.tsv", vocab_tf, vocab_gene)
            null[i] = top_keys(OUT / f"null_{i}" / "tf_to_gene.tsv", vocab_tf, vocab_gene)
            print(f"ranked run {i}  tfs {len(vocab_tf)} genes {len(vocab_gene)}", flush=True)
        np.savez(cache, real=real, null=null)

    rows = []
    hists = np.zeros((len(KS), 2, N_RUNS + 1), dtype=np.int64)
    for j, k in enumerate(KS):
        rh = support_hist(real, int(k))
        nh = support_hist(null, int(k))
        hists[j, 0] = rh
        hists[j, 1] = nh
        m, fdr, n_real, n_null = first_fdr(rh, nh)
        # also the twice-or-more bar used at the eRegulon point
        real_ge2 = int(rh[2:].sum())
        null_ge2 = int(nh[2:].sum())
        rows.append(
            {
                "k": int(k),
                "fraction_of_complete": int(k) / N_COMPLETE,
                "real_links_seen": int(rh[1:].sum()),
                "null_links_seen": int(nh[1:].sum()),
                "min_runs_fdr05": m if m is not None else np.nan,
                "fdr": fdr if fdr is not None else np.nan,
                "real_links_kept": n_real if n_real is not None else 0,
                "null_links_kept": n_null if n_null is not None else 0,
                "real_links_in_ge2": real_ge2,
                "null_links_in_ge2": null_ge2,
                "fdr_ge2": (null_ge2 / real_ge2) if real_ge2 else np.nan,
            }
        )
        print(
            f"k={k:>8}  density={k / N_COMPLETE:.5f}  "
            f"FDR<=0.05 at >={m} runs  kept {n_real} real / {n_null} null  "
            f"fdr_ge2={rows[-1]['fdr_ge2']:.4f}",
            flush=True,
        )

    table = pd.DataFrame(rows)
    table.to_csv(DEST / "sparsity_fdr.tsv", sep="\t", index=False)
    np.save(DEST / "sparsity_support_hists.npy", hists)
    print(table.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
