#!/usr/bin/env python3
"""Paired NestBoot tests for TF–gene links and eRegulons.

Each current run is one pair: real GRNBoost2 versus the matched expression-shuffle
null. For a feature, n10 is runs where it is present only in real, n01 only in
null. The exact McNemar test is a one-sided binomial test on the discordant
pairs. Benjamini–Hochberg FDR is applied within each feature class.
"""
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom, false_discovery_control, fisher_exact

LEGACY_BEFORE = datetime(2026, 10, 19, 18, 0).timestamp()
OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
N_RUNS = 50


def current_runs():
    """Runs whose real and null eRegulon tables both exist. Legacy copies are not read."""
    runs = []
    for i in range(N_RUNS):
        real = OUT / f"real_{i}" / "eRegulons.tsv"
        null = OUT / f"null_{i}" / "eRegulons.tsv"
        if real.is_file() and real.stat().st_size > 0 and null.is_file() and null.stat().st_size > 0:
            runs.append(i)
    return runs


def load_sets(kind, i):
    df = pd.read_csv(
        OUT / f"{kind}_{i}" / "eRegulons.tsv",
        sep="\t",
        usecols=["TF", "Gene", "eRegulon_name"],
    )
    links = set(zip(df["TF"].to_numpy(), df["Gene"].to_numpy()))
    ereg_genes = df.groupby("eRegulon_name")["Gene"].nunique().to_dict()
    return links, ereg_genes


def mcnemar_greater(n10, n01):
    """One-sided exact McNemar p-values. Vector inputs."""
    n10 = np.asarray(n10, dtype=np.int32)
    n01 = np.asarray(n01, dtype=np.int32)
    disc = n10 + n01
    p = np.ones(len(n10), dtype=np.float64)
    ok = disc > 0
    # P(X >= n10) under Bin(disc, 0.5)
    p[ok] = binom.sf(n10[ok] - 1, disc[ok], 0.5)
    return p


def bh(p):
    p = np.asarray(p, dtype=np.float64)
    q = np.ones_like(p)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = false_discovery_control(p[ok], method="bh")
    return q


def fisher_greater_rows(real_n, null_n, n):
    """One-sided Fisher exact p, real inclusion greater than null. Slow; call sparingly."""
    out = np.empty(len(real_n), dtype=np.float64)
    for i, (a, b) in enumerate(zip(real_n, null_n)):
        # rows: real, null; cols: present, absent
        out[i] = fisher_exact([[int(a), n - int(a)], [int(b), n - int(b)]], alternative="greater")[1]
    return out


def build_counts(runs):
    link = {}
    ereg = {}
    ereg_gene_sum = {}

    def touch(store, key, real_hit, null_hit):
        rec = store.get(key)
        if rec is None:
            rec = [0, 0, 0]  # real, null, both
            store[key] = rec
        if real_hit:
            rec[0] += 1
        if null_hit:
            rec[1] += 1
        if real_hit and null_hit:
            rec[2] += 1

    for i in runs:
        real_links, real_ereg = load_sets("real", i)
        null_links, null_ereg = load_sets("null", i)
        for key in real_links | null_links:
            touch(link, key, key in real_links, key in null_links)
        for name in set(real_ereg) | set(null_ereg):
            touch(ereg, name, name in real_ereg, name in null_ereg)
        for name, n_genes in real_ereg.items():
            slot = ereg_gene_sum.setdefault(name, [0, 0, 0, 0])  # real_sum, real_n, null_sum, null_n
            slot[0] += n_genes
            slot[1] += 1
        for name, n_genes in null_ereg.items():
            slot = ereg_gene_sum.setdefault(name, [0, 0, 0, 0])
            slot[2] += n_genes
            slot[3] += 1
        print(f"run {i}: link keys {len(link)} ereg keys {len(ereg)}", flush=True)
    return link, ereg, ereg_gene_sum


def table_from_counts(counts, n, key_names):
    rows = []
    for key, (real_n, null_n, both_n) in counts.items():
        n10 = real_n - both_n
        n01 = null_n - both_n
        rec = dict(zip(key_names, key if isinstance(key, tuple) else (key,)))
        rec.update(
            {
                "real_n": real_n,
                "null_n": null_n,
                "both_n": both_n,
                "real_only": n10,
                "null_only": n01,
                "real_frac": real_n / n,
                "null_frac": null_n / n,
            }
        )
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["p_mcnemar"] = mcnemar_greater(df["real_only"], df["null_only"])
    df["q_mcnemar"] = bh(df["p_mcnemar"].to_numpy())
    return df


def main():
    DEST.mkdir(exist_ok=True)
    runs = current_runs()
    n = len(runs)
    print(f"current runs ({n}): {runs}", flush=True)
    link_counts, ereg_counts, ereg_genes = build_counts(runs)
    links = table_from_counts(link_counts, n, ["TF", "Gene"])
    eregs = table_from_counts(ereg_counts, n, ["eRegulon"])

    gene_rows = []
    for name, (rs, rn, ns, nn) in ereg_genes.items():
        gene_rows.append(
            {
                "eRegulon": name,
                "mean_genes_real": (rs / rn) if rn else np.nan,
                "mean_genes_null": (ns / nn) if nn else np.nan,
            }
        )
    eregs = eregs.merge(pd.DataFrame(gene_rows), on="eRegulon", how="left")

    # Fisher as a pairing-free check on eRegulons (hundreds, not tens of thousands).
    eregs["p_fisher"] = fisher_greater_rows(eregs["real_n"], eregs["null_n"], n)
    eregs["q_fisher"] = bh(eregs["p_fisher"].to_numpy())

    links = links.sort_values(["q_mcnemar", "real_n", "TF", "Gene"], ascending=[True, False, True, True])
    eregs = eregs.sort_values(["q_mcnemar", "real_n", "eRegulon"], ascending=[True, False, True])

    links.to_csv(DEST / "tf_gene_nestboot.tsv", sep="\t", index=False)
    eregs.to_csv(DEST / "eregulon_nestboot.tsv", sep="\t", index=False)

    def call_counts(df, qcol="q_mcnemar"):
        sig = df[df[qcol] <= 0.05]
        return {
            "tested": len(df),
            "q05": len(sig),
            "q05_real_gt": int((sig["real_n"] > sig["null_n"]).sum()),
            "q01": int((df[qcol] <= 0.01).sum()),
        }

    link_c = call_counts(links)
    ereg_c = call_counts(eregs)
    ereg_f = call_counts(eregs, "q_fisher")
    sig_links = links[links["q_mcnemar"] <= 0.05]
    # Stability among discoveries: present in at least half of real runs.
    half = int(np.ceil(n / 2))
    stable = sig_links[sig_links["real_n"] >= half]
    sig_ereg = eregs[eregs["q_mcnemar"] <= 0.05].sort_values(
        ["real_frac", "mean_genes_real"], ascending=[False, False]
    )

    tf_counts = (
        stable.groupby("TF")
        .size()
        .sort_values(ascending=False)
        .rename("n_genes")
        .reset_index()
    )
    tf_counts.to_csv(DEST / "tf_significant_targets.tsv", sep="\t", index=False)

    # Per-run discovery is not defined; report global size test reference via inclusion.
    lines = []
    lines.append("# NestBoot tests")
    lines.append("")
    lines.append(f"Tests use all {n} pairs (runs {runs[0]}–{runs[-1]}). Every null `eRegulons.tsv` is from the current eGRN pass; 19 early nulls were rebuilt from their saved `tf_to_gene.tsv` on 22 Oct 2026. The older tables remain beside them as `eRegulons.legacy.tsv` and are not read.")
    lines.append("")
    lines.append("A feature is scored once per run, present or absent. Run i real is paired with run i null. The exact McNemar test asks whether, among runs where the two disagree, the feature is more often present in real. p-values are one-sided. Benjamini–Hochberg FDR is computed separately for TF–gene links and for eRegulon names. Absence from every tested run is not a test.")
    lines.append("")
    lines.append("## TF–gene links")
    lines.append("")
    lines.append(f"Tested **{link_c['tested']:,}** links that appear in at least one current real or null run. **{link_c['q05']:,}** have McNemar q ≤ 0.05 ({link_c['q01']:,} at q ≤ 0.01). All {link_c['q05_real_gt']:,} of the q ≤ 0.05 calls are in the real-greater direction.")
    lines.append("")
    lines.append(f"Of the q ≤ 0.05 links, **{len(stable):,}** are also present in at least {half}/{n} real runs. Those are the links counted per TF in `tf_significant_targets.tsv` ({tf_counts['TF'].nunique() if len(tf_counts) else 0} TFs). Median real support among these stable discoveries is {stable['real_n'].median() if len(stable) else float('nan'):.0f}/{n}; median null support is {stable['null_n'].median() if len(stable) else float('nan'):.0f}/{n}.")
    lines.append("")
    lines.append("| TF | significant target genes (q ≤ 0.05 and real support ≥ half) |")
    lines.append("|---|---:|")
    for row in tf_counts.head(15).itertuples(index=False):
        lines.append(f"| {row.TF} | {row.n_genes} |")
    lines.append("")
    lines.append("## eRegulons")
    lines.append("")
    lines.append(f"Tested **{ereg_c['tested']:,}** eRegulon names. **{ereg_c['q05']:,}** have McNemar q ≤ 0.05 ({ereg_c['q01']:,} at q ≤ 0.01). The unpaired Fisher exact test, same direction, calls **{ereg_f['q05']:,}** at q ≤ 0.05. Agreement on the McNemar discoveries: {(sig_ereg['q_fisher'] <= 0.05).sum() if len(sig_ereg) else 0} of {len(sig_ereg)} also have Fisher q ≤ 0.05.")
    lines.append("")
    lines.append("eRegulons with McNemar q ≤ 0.05, ordered by real inclusion then mean target-gene count:")
    lines.append("")
    lines.append("| eRegulon | real runs | null runs | mean genes real | mean genes null | q McNemar | q Fisher |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in sig_ereg.itertuples(index=False):
        ngn = "" if pd.isna(row.mean_genes_null) else f"{row.mean_genes_null:.1f}"
        ngr = "" if pd.isna(row.mean_genes_real) else f"{row.mean_genes_real:.1f}"
        lines.append(
            f"| {row.eRegulon} | {row.real_n} | {row.null_n} | {ngr} | {ngn} | {row.q_mcnemar:.2e} | {row.q_fisher:.2e} |"
        )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `results_summary/tf_gene_nestboot.tsv` — every tested TF–gene link, with discordant counts, McNemar p and q.")
    lines.append("- `results_summary/eregulon_nestboot.tsv` — every tested eRegulon, McNemar and Fisher.")
    lines.append("- `results_summary/tf_significant_targets.tsv` — per-TF count of stable real-enriched target genes.")
    lines.append("")
    lines.append("Inclusion ignores repeated regions for the same TF–gene pair inside one run. Target-gene means for an eRegulon are averaged only over runs in which that eRegulon is present.")
    text = "\n".join(lines) + "\n"
    (DEST / "NESTBOOT.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
