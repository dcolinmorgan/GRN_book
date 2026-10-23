#!/usr/bin/env python3
"""NestBoot on the raw SCENIC+ TF–gene tables, plus fixed degree cuts.

Raw: every edge written to tf_to_gene.tsv, no importance cutoff.
Degree: within each run, keep the 3, 4, or 5 strongest edges incident to
each gene (in-degree) or each TF (out-degree).

FDR at support m is
    #{shuffled links in at least m runs} / #{measured links in at least m runs}.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
PLOTS = DEST / "plots"
N_RUNS = 50
N_COMPLETE = 824 * 36601
DEGREES = (3, 4, 5)


def encode(values, vocab):
    out = np.empty(len(values), dtype=np.int32)
    for i, v in enumerate(values):
        code = vocab.get(v)
        if code is None:
            code = len(vocab)
            vocab[v] = code
        out[i] = code
    return out


def load_tf2g(path, vocab_tf, vocab_gene):
    df = pd.read_csv(path, sep="\t", usecols=["TF", "target", "importance"])
    tf = encode(df["TF"].to_numpy(), vocab_tf)
    gene = encode(df["target"].to_numpy(), vocab_gene)
    imp = df["importance"].to_numpy()
    return tf, gene, imp


def top_per_group(group, other, imp, n):
    """Strongest n edges for each group id. Returns (group, other) pairs."""
    order = np.lexsort((other, -imp, group))
    g = group[order]
    o = other[order]
    change = np.flatnonzero(g[1:] != g[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(g)]
    take = np.concatenate([np.arange(s, min(s + n, e)) for s, e in zip(starts, ends)])
    return g[take], o[take]


def _group_cuts(group):
    change = np.flatnonzero(group[1:] != group[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(group)]
    return starts, ends


def _take_first(cuts, d):
    starts, ends = cuts
    mask = np.zeros(int(ends[-1]) if len(ends) else 0, dtype=bool)
    for s, e in zip(starts, ends):
        mask[s : min(s + d, e)] = True
    return mask


def pack(tf, gene):
    return tf.astype(np.int64) * 100_000 + gene.astype(np.int64)


def support_hist(arrays):
    arrays = [np.unique(a) for a in arrays]
    stacked = np.concatenate(arrays)
    stacked.sort()
    if len(stacked) == 0:
        return np.zeros(N_RUNS + 1, dtype=np.int64)
    change = np.flatnonzero(stacked[1:] != stacked[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(stacked)]
    support = ends - starts
    return np.bincount(support, minlength=N_RUNS + 1).astype(np.int64)


def tail_table(real_hist, null_hist):
    real_ge = np.cumsum(real_hist[::-1])[::-1]
    null_ge = np.cumsum(null_hist[::-1])[::-1]
    rows = []
    for m in range(1, N_RUNS + 1):
        rr, nn = int(real_ge[m]), int(null_ge[m])
        rows.append(
            {
                "m": m,
                "measured_at": int(real_hist[m]),
                "shuffled_at": int(null_hist[m]),
                "measured_ge": rr,
                "shuffled_ge": nn,
                "fdr": (nn / rr) if rr else np.nan,
            }
        )
    return rows


def first_m(rows):
    for row in rows:
        if row["measured_ge"] > 0 and row["fdr"] <= 0.05:
            return row
    return None


def write_block(fh, title, note, edge_counts, rows):
    fh.write("\n" + "=" * 88 + "\n")
    fh.write(title + "\n")
    fh.write(note + "\n")
    fh.write(
        f"edges per run: min {min(edge_counts)}  median {int(np.median(edge_counts))}  "
        f"max {max(edge_counts)}\n"
    )
    hit = first_m(rows)
    if hit is None:
        fh.write("No support threshold reaches FDR <= 0.05.\n")
    else:
        fh.write(
            f"First m with FDR <= 0.05: links in >= {hit['m']}/50 runs  "
            f"measured {hit['measured_ge']}  shuffled {hit['shuffled_ge']}  "
            f"FDR {hit['fdr']:.6f}\n"
        )
    fh.write(
        f"{'m':>4}  {'measured_at':>12}  {'shuffled_at':>12}  "
        f"{'measured_ge':>12}  {'shuffled_ge':>12}  {'fdr_ge':>10}\n"
    )
    for row in rows:
        fdr = "nan" if np.isnan(row["fdr"]) else f"{row['fdr']:.6f}"
        fh.write(
            f"{row['m']:>4}  {row['measured_at']:>12}  {row['shuffled_at']:>12}  "
            f"{row['measured_ge']:>12}  {row['shuffled_ge']:>12}  {fdr:>10}\n"
        )


def plot_profiles(name, title, series):
    """series: list of (label, rows)."""
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 200,
        }
    )
    colors = ["#1f4e79", "#c47b2b", "#2a9d8f", "#9b2226"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))
    for (label, rows), color in zip(series, colors):
        m = [r["m"] for r in rows]
        axes[0].plot(m, [r["measured_at"] for r in rows], color=color, lw=1.6, label=label)
        axes[0].plot(
            m, [r["shuffled_at"] for r in rows], color=color, lw=1.2, ls="--", alpha=0.85
        )
        axes[1].plot(m, [r["fdr"] for r in rows], color=color, lw=1.6, label=label)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Runs in which the link is present")
    axes[0].set_ylabel("Links at this support")
    axes[0].set_title(title)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(True, ls="--", alpha=0.3)
    axes[1].axhline(0.05, color="gray", ls=":", lw=1)
    axes[1].set_xlabel("Minimum runs (m)")
    axes[1].set_ylabel("Shuffled ≥ m  /  measured ≥ m")
    axes[1].set_ylim(0, 1.5)
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(True, ls="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / f"{name}.png")
    fig.savefig(PLOTS / f"{name}.pdf")
    plt.close(fig)


def eregulon_keys():
    real, null = [], []
    counts_r, counts_n = [], []
    vocab_tf, vocab_gene = {}, {}
    for i in range(N_RUNS):
        for kind, bucket, counts in (
            ("real", real, counts_r),
            ("null", null, counts_n),
        ):
            df = pd.read_csv(
                OUT / f"{kind}_{i}" / "eRegulons.tsv",
                sep="\t",
                usecols=["TF", "Gene"],
            )
            df = df.drop_duplicates()
            tf = encode(df["TF"].to_numpy(), vocab_tf)
            gene = encode(df["Gene"].to_numpy(), vocab_gene)
            keys = pack(tf, gene)
            bucket.append(keys)
            counts.append(len(keys))
        print(f"eRegulon run {i}", flush=True)
    return real, null, counts_r, counts_n


def main():
    DEST.mkdir(exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    vocab_tf, vocab_gene = {}, {}
    raw_real, raw_null = [], []
    raw_n_real, raw_n_null = [], []
    # degree -> list over runs of key arrays
    gene_real = {d: [] for d in DEGREES}
    gene_null = {d: [] for d in DEGREES}
    tf_real = {d: [] for d in DEGREES}
    tf_null = {d: [] for d in DEGREES}
    gene_counts = {("real", d): [] for d in DEGREES}
    gene_counts.update({("null", d): [] for d in DEGREES})
    tf_counts = {("real", d): [] for d in DEGREES}
    tf_counts.update({("null", d): [] for d in DEGREES})

    for i in range(N_RUNS):
        for kind, raw_bucket, raw_counts, gene_bucket, tf_bucket in (
            ("real", raw_real, raw_n_real, gene_real, tf_real),
            ("null", raw_null, raw_n_null, gene_null, tf_null),
        ):
            tf, gene, imp = load_tf2g(OUT / f"{kind}_{i}" / "tf_to_gene.tsv", vocab_tf, vocab_gene)
            raw_bucket.append(pack(tf, gene))
            raw_counts.append(len(tf))
            # Ranked strongest-first within each gene, then within each TF.
            # The first d of a group are its top d, so one ranking serves 3, 4, and 5.
            gene_g, gene_tf = top_per_group(gene, tf, imp, max(DEGREES))
            tf_tf, tf_gene = top_per_group(tf, gene, imp, max(DEGREES))
            gene_keys = pack(gene_tf, gene_g)
            tf_keys = pack(tf_tf, tf_gene)
            gene_cut = _group_cuts(gene_g)
            tf_cut = _group_cuts(tf_tf)
            for d in DEGREES:
                gsel = _take_first(gene_cut, d)
                gene_bucket[d].append(gene_keys[gsel])
                gene_counts[(kind, d)].append(int(gsel.sum()))
                tsel = _take_first(tf_cut, d)
                tf_bucket[d].append(tf_keys[tsel])
                tf_counts[(kind, d)].append(int(tsel.sum()))
        print(
            f"tf_to_gene run {i}  edges real {raw_n_real[-1]} null {raw_n_null[-1]}  "
            f"tfs {len(vocab_tf)} genes {len(vocab_gene)}",
            flush=True,
        )

    blocks = []

    def add(title, note, counts, real_arrays, null_arrays):
        rh = support_hist(real_arrays)
        nh = support_hist(null_arrays)
        rows = tail_table(rh, nh)
        blocks.append((title, note, counts, rows))
        hit = first_m(rows)
        print(title, "first", None if hit is None else (hit["m"], hit["measured_ge"], hit["fdr"]), flush=True)
        return rows

    raw_rows = add(
        "RAW tf_to_gene: every edge SCENIC+ wrote, no importance cutoff",
        "An edge counts as present in a run when it appears in that run's tf_to_gene.tsv. "
        f"Complete bipartite graph is {N_COMPLETE:,} TF–gene pairs "
        f"({len(vocab_tf)} TFs scored, {len(vocab_gene)} genes scored).",
        raw_n_real,
        raw_real,
        raw_null,
    )
    gene_rows = {}
    tf_rows = {}
    for d in DEGREES:
        gene_rows[d] = add(
            f"SPARSITY: {d} strongest TFs per gene",
            "Within each run, each gene keeps its d highest-importance TF links. "
            "Genes with fewer than d scored edges keep all of them.",
            gene_counts[("real", d)],
            gene_real[d],
            gene_null[d],
        )
        tf_rows[d] = add(
            f"SPARSITY: {d} strongest genes per TF",
            "Within each run, each TF keeps its d highest-importance target genes. "
            "TFs with fewer than d scored edges keep all of them.",
            tf_counts[("real", d)],
            tf_real[d],
            tf_null[d],
        )

    er, en, ecr, _ecn = eregulon_keys()
    e_rows = add(
        "SCENIC+ eRegulon TF–gene links, as written (no extra cutoff)",
        "Presence in eRegulons.tsv. This is the network used in the earlier recurrence figures. "
        "It is already much smaller than the raw GRNBoost2 table.",
        ecr,
        er,
        en,
    )

    path = DEST / "scenicplus_support_fdr.txt"
    with path.open("w") as fh:
        fh.write("NestBoot support tables for the SCENIC+ benchmark\n")
        fh.write("50 measured runs and 50 gene-shuffled runs.\n")
        fh.write(
            "m is the minimum number of runs. measured_at / shuffled_at are links with "
            "support exactly m. measured_ge / shuffled_ge are links with support >= m. "
            "fdr_ge = shuffled_ge / measured_ge.\n"
        )
        fh.write(
            "The reported cutoff is the smallest m with fdr_ge <= 0.05. "
            "The ratio is not monotone, so the full column is the result.\n"
        )
        for title, note, counts, rows in blocks:
            write_block(fh, title, note, counts, rows)
    print("wrote", path, flush=True)

    plot_profiles(
        "raw_tf_to_gene_support",
        "Raw TF–gene output (solid measured, dashed shuffled)",
        [("raw tf_to_gene", raw_rows)],
    )
    plot_profiles(
        "links_per_gene",
        "3–5 TFs per gene (solid measured, dashed shuffled)",
        [(f"{d} TFs/gene", gene_rows[d]) for d in DEGREES],
    )
    plot_profiles(
        "links_per_tf",
        "3–5 genes per TF (solid measured, dashed shuffled)",
        [(f"{d} genes/TF", tf_rows[d]) for d in DEGREES],
    )
    plot_profiles(
        "eregulon_and_raw_fdr",
        "FDR curves",
        [
            ("raw tf_to_gene", raw_rows),
            ("3 TFs/gene", gene_rows[3]),
            ("5 TFs/gene", gene_rows[5]),
            ("eRegulon", e_rows),
        ],
    )


if __name__ == "__main__":
    main()
