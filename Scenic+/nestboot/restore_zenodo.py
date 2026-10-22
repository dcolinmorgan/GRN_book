#!/usr/bin/env python3
"""Restore Zenodo's flat filenames to the paths nestboot.py expects.

Zenodo bucket keys cannot contain directories, so upload_zenodo.py stores
`results_full_scenicplus/real_0/eRegulons.tsv` as
`results_full_scenicplus__real_0__eRegulons.tsv`.

Download every file from a deposit into one directory, put
zenodo_manifest.tsv next to them (or pass --manifest), then:

    python restore_zenodo.py --src downloads/inputs --dest . --deposit inputs
    python restore_zenodo.py --src downloads/outputs --dest . --deposit outputs
"""
import argparse
import hashlib
import shutil
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", required=True, help="Directory of downloaded Zenodo files")
    parser.add_argument("--dest", required=True, help="Benchmark directory to populate")
    parser.add_argument("--manifest", default=None, help="Defaults to <src>/zenodo_manifest.tsv")
    parser.add_argument("--deposit", choices=["inputs", "outputs", "both"], default="both")
    parser.add_argument("--check", action="store_true", help="Verify sha256 while copying")
    args = parser.parse_args()

    src = Path(args.src)
    dest = Path(args.dest)
    manifest = Path(args.manifest) if args.manifest else src / "zenodo_manifest.tsv"
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")

    wanted = {"inputs", "outputs"} if args.deposit == "both" else {args.deposit}
    rows = []
    for line in manifest.read_text().splitlines():
        if not line or line.startswith("deposit\t"):
            continue
        deposit, remote, rel, nbytes, digest = line.split("\t")
        if deposit in wanted:
            rows.append((deposit, remote, rel, int(nbytes), digest))
    if not rows:
        raise SystemExit(f"no rows for deposit={args.deposit}")

    missing = []
    for deposit, remote, rel, nbytes, digest in rows:
        src_file = src / remote
        if not src_file.is_file():
            missing.append(remote)
            continue
        if src_file.stat().st_size != nbytes:
            raise SystemExit(f"size mismatch for {remote}: {src_file.stat().st_size} != {nbytes}")
        if args.check and sha256(src_file) != digest:
            raise SystemExit(f"sha256 mismatch for {remote}")
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.resolve() != src_file.resolve():
            shutil.copy2(src_file, target)
        print(rel, flush=True)
    if missing:
        print("missing downloads:", file=sys.stderr)
        for name in missing:
            print(" ", name, file=sys.stderr)
        sys.exit(1)
    print(f"restored {len(rows)} files into {dest}")


if __name__ == "__main__":
    main()
