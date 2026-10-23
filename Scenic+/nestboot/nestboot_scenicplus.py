#!/usr/bin/env python3
"""NestBoot SCENIC+ benchmark, from prepared inputs through the statistics.

One real/null pair:
  real: GRNBoost2 TF_to_gene on the observed multiome (seed = --seed + run)
  null: the same after shuffling each gene's expression independently
        (seed = --seed + run + 10000)
  both: eGRN with the same region_to_gene table, direct cistromes, and
        cisTarget ranking database

region_to_gene is deterministic given the multiome and the SCENIC+ search
space, so it is an input, not part of the bootstrap. This repository does
not contain that search space; pass the table that was used for the
published runs (real_0/region_to_gene.tsv).

A side of a run is skipped when both tf_to_gene.tsv and eRegulons.tsv are
already present and non-empty. If only tf_to_gene.tsv exists, only eGRN is
repeated.
"""
import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy.sparse as sp

SCENICPLUS_COMMIT = "840dab85de3044846234c157e327b6ea0290abfb"


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def expression_modality(mudata_path: Path) -> str:
    import mudata

    mdata = mudata.read(mudata_path)
    for key in ("scRNA", "rna", "GEX"):
        if key in mdata.mod:
            return key
    raise SystemExit(f"No expression modality in {mudata_path}; looked for scRNA, rna, GEX")


def create_shuffled_mudata(mdata_path, expr_mod, out_path, seed):
    """MuData with each gene's expression shuffled independently."""
    import anndata
    import mudata

    np.random.seed(seed)
    mdata = mudata.read(mdata_path)
    expr = mdata.mod[expr_mod]
    X = expr.X.toarray() if sp.issparse(expr.X) else np.array(expr.X)
    for j in range(X.shape[1]):
        np.random.shuffle(X[:, j])
    shuffled = anndata.AnnData(X=sp.csr_matrix(X), obs=expr.obs.copy(), var=expr.var.copy())
    mudata.MuData({expr_mod: shuffled, **{k: v for k, v in mdata.mod.items() if k != expr_mod}}).write(out_path)
    return out_path


def scenicplus(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(result.stderr[-3000:])
        return False
    return True


def run_tf_to_gene(mudata_path, tf_names, out_dir, n_cpu, seed):
    out_dir = Path(out_dir)
    tmp = out_dir / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "tf_to_gene.tsv"
    ok = scenicplus([
        "scenicplus", "grn_inference", "TF_to_gene",
        "--multiome_mudata_fname", str(mudata_path),
        "--tf_names", str(tf_names),
        "--temp_dir", str(tmp),
        "--out_tf_to_gene_adjacencies", str(dest),
        "--method", "GBM",
        "--n_cpu", str(n_cpu),
        "--seed", str(seed),
    ])
    return dest if ok else None


def run_egrn(tf2g, region_to_gene, cistromes, ranking_db, out_dir, n_cpu, seed):
    out_dir = Path(out_dir)
    tmp = out_dir / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "eRegulons.tsv"
    ok = scenicplus([
        "scenicplus", "grn_inference", "eGRN",
        "--TF_to_gene_adj_fname", str(tf2g),
        "--region_to_gene_adj_fname", str(region_to_gene),
        "--cistromes_fname", str(cistromes),
        "--ranking_db_fname", str(ranking_db),
        "--eRegulon_out_fname", str(dest),
        "--temp_dir", str(tmp),
        "--n_cpu", str(n_cpu),
        "--seed", str(seed),
    ])
    return dest if ok else None


def run_side(kind, run, mudata_path, tf_names, region_to_gene, cistromes, ranking_db, out_dir, n_cpu, seed, expr_mod, is_null):
    dest = Path(out_dir) / f"{kind}_{run}"
    dest.mkdir(parents=True, exist_ok=True)
    tf2g = dest / "tf_to_gene.tsv"
    ereg = dest / "eRegulons.tsv"
    if nonempty(tf2g) and nonempty(ereg):
        log(f"  {kind} run {run}: already complete, skipping")
        return
    if not nonempty(tf2g):
        source = mudata_path
        if is_null:
            log(f"  NULL run {run}: shuffling expression (seed {seed})")
            source = dest / "shuffled.h5mu"
            create_shuffled_mudata(mudata_path, expr_mod, source, seed)
        log(f"  {kind} run {run}: TF_to_gene (seed {seed}, n_cpu {n_cpu})")
        written = run_tf_to_gene(source, tf_names, dest, n_cpu, seed)
        if is_null and Path(source).exists():
            Path(source).unlink()
        if written is None:
            raise SystemExit(f"TF_to_gene failed for {kind}_{run}")
        tf2g = Path(written)
    else:
        log(f"  {kind} run {run}: reusing {tf2g.name}")
    log(f"  {kind} run {run}: eGRN (seed {seed})")
    written = run_egrn(tf2g, region_to_gene, cistromes, ranking_db, dest, n_cpu, seed)
    if written is None or not nonempty(Path(written)):
        raise SystemExit(f"eGRN failed for {kind}_{run}")
    log(f"  {kind} run {run}: wrote {written}")


def cmd_run(args):
    for label, path in (
        ("mudata", args.mudata),
        ("tf_names", args.tf_names),
        ("region_to_gene", args.region_to_gene),
        ("cistromes", args.cistromes),
        ("ranking_db", args.ranking_db),
    ):
        if not Path(path).is_file():
            raise SystemExit(f"Missing {label}: {path}")
    expr_mod = expression_modality(Path(args.mudata))
    log(f"scenicplus pin: {SCENICPLUS_COMMIT}")
    log(f"Expression modality: {expr_mod}")
    log(f"Runs {args.start_run}..{args.end_run - 1}, n_cpu={args.n_cpu}, base seed={args.seed}")
    for i in range(args.start_run, args.end_run):
        log("=" * 60)
        log(f"Run {i} ({i - args.start_run + 1}/{args.end_run - args.start_run})")
        run_side(
            "real", i, args.mudata, args.tf_names, args.region_to_gene,
            args.cistromes, args.ranking_db, args.out_dir, args.n_cpu,
            args.seed + i, expr_mod, False,
        )
        run_side(
            "null", i, args.mudata, args.tf_names, args.region_to_gene,
            args.cistromes, args.ranking_db, args.out_dir, args.n_cpu,
            args.seed + i + 10000, expr_mod, True,
        )
    log("BOOTSTRAP COMPLETE")


def cmd_stats(_args):
    import build_nestboot_stats
    import nestboot_analysis

    build_nestboot_stats.main()
    nestboot_analysis.main()


def main():
    parser = argparse.ArgumentParser(description="NestBoot SCENIC+ benchmark")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run real/null bootstraps. Skips finished sides.")
    run.add_argument("--mudata", default="work/nestboot_scenicplus/base_mudata.h5mu")
    run.add_argument("--tf_names", default="work/menr/tf_names.txt")
    run.add_argument("--region_to_gene", default="results_full_scenicplus/real_0/region_to_gene.tsv")
    run.add_argument("--cistromes", default="work/menr/cistromes_direct.h5ad")
    run.add_argument("--ranking_db", default="hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")
    run.add_argument("--out_dir", default="results_full_scenicplus")
    run.add_argument("--start_run", type=int, default=0)
    run.add_argument("--end_run", type=int, default=50)
    run.add_argument("--n_cpu", type=int, default=8)
    run.add_argument("--seed", type=int, default=42)
    run.set_defaults(func=cmd_run)

    stats = sub.add_parser("stats", help="Rebuild eRegulon size summaries from eRegulons.tsv")
    stats.set_defaults(func=cmd_stats)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
