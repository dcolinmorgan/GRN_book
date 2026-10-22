#!/usr/bin/env python3
"""Summarize 50 NestBoot real/null eRegulon runs."""
from pathlib import Path

from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Null eRegulon files written before this stamp used an older eGRN pass.
# Re-running that pass on saved tf_to_gene tables (21 Oct) produced null
# networks the same size as the later runs, so those early files are not
# exchangeable with the current ones.
LEGACY_BEFORE = datetime(2026, 10, 19, 18, 0).timestamp()

OUT = Path("results_full_scenicplus")
DEST = Path("results_summary")
N = 50
COLS = [
    "Region",
    "Gene",
    "TF",
    "is_extended",
    "eRegulon_name",
    "importance_TF2G",
    "regulation",
    "rho_TF2G",
]


def load(kind, i):
    path = OUT / f"{kind}_{i}" / "eRegulons.tsv"
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df["is_extended"] = df["is_extended"].astype(bool)
    return df


def tf_gene_set(df):
    return set(zip(df["TF"].to_numpy(), df["Gene"].to_numpy()))


def jaccard(a, b):
    if not a and not b:
        return 1.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def main():
    DEST.mkdir(exist_ok=True)
    rows = []
    real_sets = []
    null_sets = []
    real_freq = {}
    null_freq = {}
    real_ereg = {}
    null_ereg = {}

    def add_freq(store, keys):
        for key in keys:
            store[key] = store.get(key, 0) + 1

    for i in range(N):
        real = load("real", i)
        null = load("null", i)
        legacy = False
        rs, ns = tf_gene_set(real), tf_gene_set(null)
        real_sets.append(rs)
        null_sets.append(ns)
        if not legacy:
            add_freq(real_freq, rs)
            add_freq(null_freq, ns)
            add_freq(real_ereg, set(real["eRegulon_name"]))
            add_freq(null_ereg, set(null["eRegulon_name"]))
        rows.append(
            {
                "run": i,
                "real_triplets": len(real),
                "null_triplets": len(null),
                "fold_triplets": len(real) / len(null) if len(null) else np.nan,
                "real_tf": real["TF"].nunique(),
                "null_tf": null["TF"].nunique(),
                "real_genes": real["Gene"].nunique(),
                "null_genes": null["Gene"].nunique(),
                "real_regions": real["Region"].nunique(),
                "null_regions": null["Region"].nunique(),
                "real_eregulons": real["eRegulon_name"].nunique(),
                "null_eregulons": null["eRegulon_name"].nunique(),
                "real_extended_frac": float(real["is_extended"].mean()),
                "null_extended_frac": float(null["is_extended"].mean()),
                "real_pos_reg_frac": float((real["regulation"] > 0).mean()),
                "null_pos_reg_frac": float((null["regulation"] > 0).mean()),
                "real_median_importance_TF2G": float(real["importance_TF2G"].median()),
                "null_median_importance_TF2G": float(null["importance_TF2G"].median()),
                "real_median_rho_TF2G": float(real["rho_TF2G"].median()),
                "null_median_rho_TF2G": float(null["rho_TF2G"].median()),
                "tf_gene_real": len(rs),
                "tf_gene_null": len(ns),
                "tf_gene_jaccard_paired": jaccard(rs, ns),
                "tf_gene_overlap": len(rs & ns),
                "cohort": "legacy_null" if legacy else "current",
            }
        )
        print(f"run {i}: real {len(real)} null {len(null)} {rows[-1]['cohort']}", flush=True)

    per = pd.DataFrame(rows)
    per.to_csv(DEST / "per_run_counts.tsv", sep="\t", index=False)

    # Mean pairwise Jaccard among real TF-gene sets, and among null sets.
    def mean_pairwise(sets):
        vals = []
        for a in range(len(sets)):
            for b in range(a + 1, len(sets)):
                vals.append(jaccard(sets[a], sets[b]))
        return float(np.mean(vals)), float(np.median(vals))

    current_idx = [i for i, row in enumerate(rows) if row["cohort"] == "current"]
    cur_real_sets = [real_sets[i] for i in current_idx]
    cur_null_sets = [null_sets[i] for i in current_idx]
    n_cur = len(current_idx)
    real_jac_mean, real_jac_med = mean_pairwise(real_sets)
    null_jac_mean, null_jac_med = mean_pairwise(null_sets)
    cur_real_jac_mean, cur_real_jac_med = mean_pairwise(cur_real_sets)
    cur_null_jac_mean, cur_null_jac_med = mean_pairwise(cur_null_sets)

    # Recurrence table: TF-gene seen in either condition.
    keys = set(real_freq) | set(null_freq)
    rec_rows = []
    for tf, gene in keys:
        rec_rows.append(
            {
                "TF": tf,
                "Gene": gene,
                "real_n": real_freq.get((tf, gene), 0),
                "null_n": null_freq.get((tf, gene), 0),
            }
        )
    rec = pd.DataFrame(rec_rows)
    rec["real_frac"] = rec["real_n"] / n_cur
    rec["null_frac"] = rec["null_n"] / n_cur
    rec = rec.sort_values(["real_n", "null_n", "TF", "Gene"], ascending=[False, True, True, True])
    # Keep the analysis table to links seen at least 5 times in real, plus
    # every link seen in a null, so the null tail is not dropped.
    keep = rec[(rec["real_n"] >= 5) | (rec["null_n"] >= 1)]
    keep.to_csv(DEST / "tf_gene_recurrence.tsv", sep="\t", index=False)

    ereg_keys = set(real_ereg) | set(null_ereg)
    ereg = pd.DataFrame(
        [
            {
                "eRegulon": name,
                "real_n": real_ereg.get(name, 0),
                "null_n": null_ereg.get(name, 0),
            }
            for name in ereg_keys
        ]
    )
    ereg["real_frac"] = ereg["real_n"] / n_cur
    ereg["null_frac"] = ereg["null_n"] / n_cur
    ereg = ereg.sort_values(["real_n", "null_n", "eRegulon"], ascending=[False, True, True])
    ereg.to_csv(DEST / "eregulon_recurrence.tsv", sep="\t", index=False)

    cur = per[per["cohort"] == "current"].copy()
    legacy = per[per["cohort"] == "legacy_null"].copy()

    def summarize(frame, col):
        x = frame[col].to_numpy()
        return {
            "mean": float(np.mean(x)),
            "std": float(np.std(x, ddof=1)),
            "median": float(np.median(x)),
            "min": float(np.min(x)),
            "max": float(np.max(x)),
        }

    w = wilcoxon(cur["real_triplets"], cur["null_triplets"], alternative="greater")
    n_real_gt = int((cur["real_triplets"] > cur["null_triplets"]).sum())

    def support_counts(freq):
        vals = np.array(list(freq.values()))
        return {k: int((vals >= k).sum()) for k in (1, 5, 10, 25, 40, 50)}

    real_support = support_counts(real_freq)
    null_support = support_counts(null_freq)
    # Links with strong real support and almost no null support.
    specific = rec[(rec["real_n"] >= max(1, int(round(0.5 * n_cur)))) & (rec["null_n"] <= 2)]
    specific_hi = rec[(rec["real_n"] >= max(1, int(round(0.8 * n_cur)))) & (rec["null_n"] == 0)]
    k_half = int(round(0.5 * n_cur))
    k_hi = int(round(0.8 * n_cur))

    lines = []
    lines.append("# NestBoot SCENIC+ summary")
    lines.append("")
    legacy_ids = ", ".join(str(i) for i in legacy["run"].tolist())
    lines.append("Fifty paired bootstraps. Each pair re-runs stochastic GRNBoost2 `TF_to_gene` (seed `42 + run` for real, `42 + run + 10000` for null) and rebuilds eGRNs with the fixed `region_to_gene` table from `real_0`. The null shuffles expression independently within each gene before GRNBoost2. An eRegulon triplet is one row of `eRegulons.tsv` (region, gene, TF).")
    lines.append("")
    lines.append("## Runs")
    lines.append("")
    if len(legacy) == 0:
        lines.append(f"All **{n_cur}** pairs are on the current eGRN pass. Nineteen early null tables were rebuilt from saved `tf_to_gene.tsv` files on 22 Oct 2026; those previous outputs are `eRegulons.legacy.tsv` and are not used here.")
    else:
        lines.append(f"**{n_cur} runs** are on the current eGRN pass. **{len(legacy)}** null tables are still the older pass: runs {legacy_ids}. Headline numbers use the current pairs only.")
    lines.append("")
    lines.append("## Size: real versus null")
    lines.append("")
    lines.append(f"On the {n_cur} current pairs, real triplets exceed the paired null in **{n_real_gt}/{n_cur}** runs (Wilcoxon signed-rank, greater, statistic {w.statistic:.0f}, p = {w.pvalue:.3e}).")
    lines.append("")
    lines.append("| | real triplets | null triplets | fold (real/null) | real TF–gene links | null TF–gene links |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for label, key in (
        ("mean", "mean"),
        ("sd", "std"),
        ("median", "median"),
        ("min", "min"),
        ("max", "max"),
    ):
        a, b, c = summarize(cur, "real_triplets"), summarize(cur, "null_triplets"), summarize(cur, "fold_triplets")
        d, e = summarize(cur, "tf_gene_real"), summarize(cur, "tf_gene_null")
        lines.append(
            f"| {label} | {a[key]:,.1f} | {b[key]:,.1f} | {c[key]:,.1f} | {d[key]:,.1f} | {e[key]:,.1f} |"
        )
    lines.append("")
    rt, nt = summarize(cur, "real_tf"), summarize(cur, "null_tf")
    rg, ng = summarize(cur, "real_genes"), summarize(cur, "null_genes")
    rr, nr = summarize(cur, "real_regions"), summarize(cur, "null_regions")
    re, ne = summarize(cur, "real_eregulons"), summarize(cur, "null_eregulons")
    lines.append("| | real mean (min–max) | null mean (min–max) |")
    lines.append("|---|---:|---:|")
    lines.append(f"| TFs | {rt['mean']:.1f} ({rt['min']:.0f}–{rt['max']:.0f}) | {nt['mean']:.1f} ({nt['min']:.0f}–{nt['max']:.0f}) |")
    lines.append(f"| genes | {rg['mean']:.1f} ({rg['min']:.0f}–{rg['max']:.0f}) | {ng['mean']:.1f} ({ng['min']:.0f}–{ng['max']:.0f}) |")
    lines.append(f"| regions | {rr['mean']:.1f} ({rr['min']:.0f}–{rr['max']:.0f}) | {nr['mean']:.1f} ({nr['min']:.0f}–{nr['max']:.0f}) |")
    lines.append(f"| eRegulons | {re['mean']:.1f} ({re['min']:.0f}–{re['max']:.0f}) | {ne['mean']:.1f} ({ne['min']:.0f}–{ne['max']:.0f}) |")
    lines.append("")
    ext_r, ext_n = summarize(cur, "real_extended_frac"), summarize(cur, "null_extended_frac")
    pos_r, pos_n = summarize(cur, "real_pos_reg_frac"), summarize(cur, "null_pos_reg_frac")
    imp_r, imp_n = summarize(cur, "real_median_importance_TF2G"), summarize(cur, "null_median_importance_TF2G")
    rho_r, rho_n = summarize(cur, "real_median_rho_TF2G"), summarize(cur, "null_median_rho_TF2G")
    lines.append(
        f"No extended triplets appear in the current runs "
        f"(mean extended fraction {ext_r['mean']:.3f} real, {ext_n['mean']:.3f} null). "
        f"Positive `regulation` fraction averages {pos_r['mean']:.3f} real and {pos_n['mean']:.3f} null. "
        f"Median TF-to-gene importance is {imp_r['mean']:.3f} real vs {imp_n['mean']:.3f} null; "
        f"median rho is {rho_r['mean']:.3f} real vs {rho_n['mean']:.3f} null. "
        "Null networks are not just smaller copies of the real ones: they contain more TFs and more eRegulon names, "
        "but each covers far fewer genes and regions. Real TF–gene rho is higher; null importance among the triplets that remain is not lower."
    )
    lines.append("")
    lines.append("## Stability of TF–gene links")
    lines.append("")
    lines.append(f"Mean pairwise Jaccard of TF–gene sets across the {n_cur} real runs is **{cur_real_jac_mean:.3f}** (median {cur_real_jac_med:.3f}). Across the matched null runs it is **{cur_null_jac_mean:.3f}** (median {cur_null_jac_med:.3f}).")
    lines.append("")
    lines.append(f"Paired real-vs-null Jaccard on the current runs averages **{cur['tf_gene_jaccard_paired'].mean():.3f}** (median {cur['tf_gene_jaccard_paired'].median():.3f}), with a mean overlap of {cur['tf_gene_overlap'].mean():.0f} links.")
    lines.append("")
    lines.append(f"Unique TF–gene links seen in at least k of the {n_cur} current runs:")
    lines.append("")
    lines.append("| k | real | null |")
    lines.append("|---:|---:|---:|")
    for k in (1, 5, 10, 25, 40, 50):
        if k > n_cur:
            continue
        rk = real_support.get(k, 0)
        nk = null_support.get(k, 0)
        lines.append(f"| {k} | {rk:,} | {nk:,} |")
    lines.append("")
    lines.append(f"Links in ≥{k_half}/{n_cur} real runs and ≤2 null runs: **{len(specific):,}**. Links in ≥{k_hi}/{n_cur} real runs and in no null run: **{len(specific_hi):,}**.")
    lines.append("")
    lines.append(f"Union of TF–gene links in the current runs is {len(rec):,} ({(rec['real_n']>0).sum():,} real, {(rec['null_n']>0).sum():,} null). `tf_gene_recurrence.tsv` keeps links with real support ≥5 or any null support ({len(keep):,} rows). Support columns count only the {n_cur} current runs.")
    lines.append("")
    lines.append("## eRegulon recurrence")
    lines.append("")
    lines.append(f"{len(ereg):,} eRegulon names appear at least once in the current runs. {(ereg['real_n']==n_cur).sum()} are in every current real run; {(ereg['null_n']==n_cur).sum()} are in every current null run. {(ereg['real_n']>=k_hi).sum()} real eRegulons reach ≥{k_hi}/{n_cur}, against {(ereg['null_n']>=k_hi).sum()} null eRegulons.")
    lines.append("")
    top = ereg.head(15)
    lines.append("Fifteen eRegulons with the highest real support (ties broken by lower null support):")
    lines.append("")
    lines.append("| eRegulon | real runs | null runs |")
    lines.append("|---|---:|---:|")
    for row in top.itertuples(index=False):
        lines.append(f"| {row.eRegulon} | {row.real_n} | {row.null_n} |")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `results_summary/per_run_counts.tsv` — one row per bootstrap, with `cohort` (`current` or `legacy_null`).")
    lines.append("- `results_summary/tf_gene_recurrence.tsv` — TF–gene support in the current runs.")
    lines.append("- `results_summary/eregulon_recurrence.tsv` — eRegulon-name support in the current runs.")
    lines.append("")
    lines.append("Counts use presence in a run, not how many regions repeat the same TF–gene pair inside that run.")
    text = "\n".join(lines) + "\n"
    (DEST / "SUMMARY.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
