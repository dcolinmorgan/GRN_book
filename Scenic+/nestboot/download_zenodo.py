#!/usr/bin/env python3
"""Download the single NestBoot Zenodo record into this directory.

The record is draft 23559287 until it is published. A draft download needs
ZENODO_TOKEN. After publication the same script uses the public record API
and the token can be omitted.

The rankings database is inside this record. There is no separate Aerts-lab
download. outputs_all.zip is unpacked to results_full_scenicplus/.
"""

import os
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen

import json

RECORD = "23559287"
DEST = Path(".")


def get(url, token):
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(url, headers=headers)) as response:
        return json.load(response)


def download(url, path, token):
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(Request(url, headers=headers)) as response, path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)


def files_from_draft(token):
    deposit = get(f"https://zenodo.org/api/deposit/depositions/{RECORD}", token)
    rows = []
    for item in deposit.get("files", []):
        rows.append((item["filename"], item["links"]["download"]))
    return rows


def files_from_record():
    record = get(f"https://zenodo.org/api/records/{RECORD}", None)
    rows = []
    for item in record.get("files", []):
        rows.append((item["key"], item["links"]["self"]))
    return rows


def main():
    token = os.environ.get("ZENODO_TOKEN", "")
    try:
        rows = files_from_record()
    except Exception:
        if not token:
            raise SystemExit("Record is not public. Set ZENODO_TOKEN to download the draft.")
        rows = files_from_draft(token)
    for name, url in rows:
        target = DEST / name
        if target.is_file() and target.stat().st_size > 0:
            print("have", name, flush=True)
            continue
        print("get", name, flush=True)
        download(url, target, token)
    archive = DEST / "outputs_all.zip"
    if archive.is_file():
        print("unpack", archive.name, flush=True)
        with zipfile.ZipFile(archive) as handle:
            handle.extractall(DEST)
    print("done", flush=True)


if __name__ == "__main__":
    main()
