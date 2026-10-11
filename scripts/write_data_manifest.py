#!/usr/bin/env python
"""
Raw-input checksum manifest — the chain-of-custody record for data that is
too large to commit (CIM 2018 §2.2 / §2.7.3: reliable storage, and a way to
prove the inputs behind a result are the ones that were archived).

Hashes every raw input the analysis reads — MRDS, ASTER TIR rasters, AVIRIS
reflectance cubes, cached fault layers — and writes one line per file:

    <sha256>  <size_bytes>  <path relative to repo root>

Derived products (zones, ratio mosaics, figures, results) are excluded: they
are regenerable from these inputs plus the git commit recorded in each
results/{site_id}_provenance.json.

Outputs
-------
    data/MANIFEST.sha256   (small text file; committed)

Usage
-----
    conda run -n aster-minerals python scripts/write_data_manifest.py
    conda run -n aster-minerals python scripts/write_data_manifest.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = REPO_ROOT / "data" / "MANIFEST.sha256"

PATTERNS = (
    "data/mrds.csv",
    "data/mrds.met",
    "data/mrds-csv.zip",
    "data/aster/*_TIR_B1[0-4].tif",
    "data/sites/*/aster/AST_L1T_*_TIR_B1[0-4].tif",
    "data/gemx/aviris/*/*.nc",
    "data/structures/*.geojson",
)


def input_files() -> list[Path]:
    files: set[Path] = set()
    for pattern in PATTERNS:
        files.update(p for p in REPO_ROOT.glob(pattern) if p.is_file())
    return sorted(files)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _entry(path: Path) -> str:
    return f"{sha256(path)}  {path.stat().st_size}  {path.relative_to(REPO_ROOT)}"


def write() -> None:
    files = input_files()
    print(f"Hashing {len(files)} files ...", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        lines = list(pool.map(_entry, files))
    MANIFEST.write_text("\n".join(lines) + "\n")
    total_gb = sum(p.stat().st_size for p in files) / 1e9
    print(f"Wrote {MANIFEST.relative_to(REPO_ROOT)} ({len(lines)} files, {total_gb:.1f} GB)")


def check() -> int:
    if not MANIFEST.is_file():
        print(f"No manifest at {MANIFEST}")
        return 2
    recorded = {}
    for line in MANIFEST.read_text().splitlines():
        digest, _size, rel = line.split("  ", 2)
        recorded[rel] = digest
    on_disk = {str(p.relative_to(REPO_ROOT)): p for p in input_files()}

    missing = sorted(set(recorded) - set(on_disk))
    new = sorted(set(on_disk) - set(recorded))
    with ThreadPoolExecutor(max_workers=4) as pool:
        common = sorted(set(recorded) & set(on_disk))
        digests = dict(zip(common, pool.map(lambda r: sha256(on_disk[r]), common)))
    changed = [r for r in common if digests[r] != recorded[r]]

    for label, items in (("MISSING", missing), ("CHANGED", changed), ("NEW (unrecorded)", new)):
        for rel in items:
            print(f"{label:17s} {rel}")
    print(f"{len(common) - len(changed)}/{len(recorded)} recorded files verified")
    return 1 if (missing or changed) else 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="verify files against the manifest")
    args = ap.parse_args()
    sys.exit(check() if args.check else write())


if __name__ == "__main__":
    main()
