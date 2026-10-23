#!/usr/bin/env python3
"""Create two Zenodo depositions: full-run inputs, and bootstrap outputs.

Requires ZENODO_TOKEN (personal access token with deposit:write).
Does not publish the records; they are left as drafts so the creator,
description, and community can be checked before "Publish".

  python upload_zenodo.py --creator "Family, Given" --dry-run
  ZENODO_TOKEN=... python upload_zenodo.py --creator "Family, Given"
"""
import argparse
import hashlib
import os
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent

INPUTS = [
    "work/nestboot_scenicplus/base_mudata.h5mu",
    "work/menr/tf_names.txt",
    "work/menr/cistromes_direct.h5ad",
    "work/menr/cistromes_extended.h5ad",
    "results_full_scenicplus/real_0/region_to_gene.tsv",
    "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
]

CODE = [
    "nestboot_scenicplus.py",
    "nestboot_optimized.py",
    "build_nestboot_stats.py",
    "nestboot_analysis.py",
    "restore_zenodo.py",
    "README.md",
    "requirements.txt",
    "requirements.lock",
    "build-constraints.txt",
]


def output_files():
    files = []
    base = ROOT / "results_full_scenicplus"
    for kind in ("real", "null"):
        for i in range(50):
            for name in ("eRegulons.tsv", "tf_to_gene.tsv"):
                files.append(f"results_full_scenicplus/{kind}_{i}/{name}")
    summary = ROOT / "results_summary"
    for path in sorted(summary.iterdir()):
        if path.is_file() and path.suffix in {".tsv", ".md"}:
            files.append(str(path.relative_to(ROOT)))
    return files


DEPOSITS = {
    "inputs": {
        "title": "NestBoot SCENIC+ benchmark inputs",
        "description": (
            "<p>Inputs required to rerun the NestBoot SCENIC+ benchmark with "
            "<code>python nestboot_scenicplus.py run</code>. scenicplus 1.0a2 commit "
            "840dab85de3044846234c157e327b6ea0290abfb. "
            "region_to_gene.tsv is the fixed region-to-gene table used for every "
            "bootstrap (it is not recomputed inside the bootstrap). "
            "The rankings file is hg38 screen v10 cluster regions vs motifs "
            "(cisTarget).</p>"
        ),
        "files": INPUTS + CODE,
    },
    "outputs": {
        "title": "NestBoot SCENIC+ benchmark outputs (50 real/null pairs)",
        "description": (
            "<p>Outputs of 50 paired NestBoot SCENIC+ runs. Each "
            "<code>real_i</code> and <code>null_i</code> directory contains "
            "TF-to-gene adjacencies and eRegulons. "
            "<code>results_summary/</code> is the NestBoot support summary, and is enough for the frequency curves without "
            "rerunning GRNBoost2. Null eRegulons were produced by the eGRN "
            "command in nestboot_scenicplus.py (scenicplus 840dab85); "
            "eRegulons.legacy.tsv files from an earlier pass are not included.</p>"
        ),
        "files": None,
    },
}


def remote_name(rel: str) -> str:
    return rel.replace("/", "__")


def load_cached_hashes(path: Path):
    cached = {}
    if not path.is_file():
        return cached
    for line in path.read_text().splitlines():
        if not line or line.startswith("deposit\t"):
            continue
        _deposit, _remote, rel, nbytes, digest = line.split("\t")
        cached[rel] = (int(nbytes), digest)
    return cached


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(selected):
    """Flat Zenodo names plus sha256. Reuses a hash when size is unchanged."""
    manifest_path = ROOT / "zenodo_manifest.tsv"
    cached = load_cached_hashes(manifest_path)
    rows = ["deposit\tzenodo_name\tpath\tbytes\tsha256"]
    for key in selected:
        for rel in DEPOSITS[key]["files"]:
            path = ROOT / rel
            size = path.stat().st_size
            previous = cached.get(rel)
            if previous and previous[0] == size:
                digest = previous[1]
            else:
                print(f"hash {rel}", flush=True)
                digest = file_sha256(path)
            rows.append(f"{key}\t{remote_name(rel)}\t{rel}\t{size}\t{digest}")
    manifest_path.write_text("\n".join(rows) + "\n")
    print(f"wrote {manifest_path}")
    return manifest_path


def check(files):
    missing = []
    total = 0
    for rel in files:
        path = ROOT / rel
        if not path.is_file():
            missing.append(rel)
            continue
        total += path.stat().st_size
    return missing, total


def put_file(bucket, token, remote, path: Path, attempts=6):
    url = f"{bucket}/{remote}"
    headers = {"Authorization": f"Bearer {token}"}
    for attempt in range(1, attempts + 1):
        try:
            with path.open("rb") as handle:
                response = requests.put(url, data=handle, headers=headers, timeout=(60, 3600))
            response.raise_for_status()
            return
        except requests.RequestException as exc:
            print(f"  attempt {attempt}/{attempts} failed for {remote}: {exc}", flush=True)
            if attempt == attempts:
                raise
            time.sleep(min(60, 5 * attempt))


def existing_files(api, token, deposition_id):
    response = requests.get(
        f"{api}/deposit/depositions/{deposition_id}",
        headers={"Authorization": f"Bearer {token}"},
        timeout=120,
    )
    response.raise_for_status()
    draft = response.json()
    present = {}
    for item in draft.get("files") or []:
        present[item.get("filename")] = item.get("filesize") or 0
    return draft, present


def create_draft(api, token, title, description, creator):
    metadata = {
        "metadata": {
            "title": title,
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": creator}],
            "keywords": ["SCENIC+", "gene regulatory network", "bootstrap", "cisTarget"],
            "license": "cc-by-4.0",
        }
    }
    response = requests.post(
        f"{api}/deposit/depositions",
        json=metadata,
        headers={"Authorization": f"Bearer {token}"},
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--creator", default="", help='Zenodo creator, "Family, Given"')
    parser.add_argument("--api", default="https://zenodo.org/api")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--manifest-only", action="store_true", help="Write zenodo_manifest.tsv and exit")
    parser.add_argument("--which", choices=["inputs", "outputs", "both"], default="both")
    parser.add_argument("--resume", type=int, default=0, help="Existing deposition id for the inputs draft")
    args = parser.parse_args()

    token = os.environ.get("ZENODO_TOKEN", "")
    DEPOSITS["outputs"]["files"] = output_files() + CODE
    selected = ["inputs", "outputs"] if args.which == "both" else [args.which]

    for key in selected:
        spec = DEPOSITS[key]
        missing, total = check(spec["files"])
        print(f"{key}: {len(spec['files'])} files, {total / 1e9:.2f} GB")
        if missing:
            print("missing:")
            for rel in missing:
                print(" ", rel)
            sys.exit(1)
    write_manifest(selected)
    if args.dry_run or args.manifest_only:
        return
    if not args.creator:
        raise SystemExit("--creator is required to create a draft")
    if not token:
        raise SystemExit("ZENODO_TOKEN is not set. Drafts were not created.")
    for key in selected:
        spec = DEPOSITS[key]
        if key == "inputs" and args.resume:
            draft, present = existing_files(args.api, token, args.resume)
        else:
            draft = create_draft(args.api, token, spec["title"], spec["description"], args.creator)
            present = {}
        bucket = draft["links"]["bucket"]
        print(f"draft {draft['id']} {draft['links']['html']}", flush=True)
        # Bucket keys are flat filenames. Prefix with the relative path so
        # real_0/eRegulons.tsv does not collide with null_0/eRegulons.tsv.
        manifest = ROOT / "zenodo_manifest.tsv"
        remote = manifest.name
        if present.get(remote) == manifest.stat().st_size:
            print(f"  skip {remote}", flush=True)
        else:
            print(f"  upload {remote}", flush=True)
            put_file(bucket, token, remote, manifest)
        for rel in spec["files"]:
            path = ROOT / rel
            remote = remote_name(rel)
            if present.get(remote) == path.stat().st_size:
                print(f"  skip {rel}", flush=True)
                continue
            print(f"  upload {rel} -> {remote}", flush=True)
            put_file(bucket, token, remote, path)
        print(f"{key} draft ready (not published): {draft['links']['html']}", flush=True)


if __name__ == "__main__":
    main()
