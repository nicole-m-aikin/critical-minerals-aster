"""Re-mosaic (with inter-granule normalisation) and regenerate figures.

Two modes
---------
--remosaic   Re-download all granules from EarthAccess, rebuild the
             normalised mosaic, then regenerate figures.  Requires a
             valid ~/.netrc entry for urs.earthdata.nasa.gov.
             Downloads are rate-limited to MAX_PARALLEL_DOWNLOADS
             concurrent sites to avoid hammering EarthAccess.

(default)    Regenerate figures only from whatever TIR data is already
             on disk (mosaic or single-granule).  Fast, no network needed.

Usage
-----
    # Re-download + re-mosaic + figures (what you want after the norm fix):
    conda run -n aster-minerals python scripts/regen_figures_parallel.py --remosaic [site_id ...]

    # Figures only (no download):
    conda run -n aster-minerals python scripts/regen_figures_parallel.py [site_id ...]

If no site IDs are given, all sites that have ASTER data on disk are processed.
"""

from __future__ import annotations

import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "sites"

# Cap concurrent EarthAccess downloads to avoid rate-limiting.
MAX_PARALLEL_DOWNLOADS = 4

# ---------------------------------------------------------------------------
# Worker scripts (run in subprocesses so each site gets its own interpreter
# and a failure in one doesn't kill the others).
# ---------------------------------------------------------------------------

_REMOSAIC_WORKER = """
import sys
from pathlib import Path

site_id  = sys.argv[1]
repo_root = Path(sys.argv[2])

import earthaccess
earthaccess.login(strategy="netrc")

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import download_and_mosaic_aster, run_site

site  = load_site_by_id(site_id, repo_root / "sites")
paths = site_paths_for(site, repo_root)

print(f"[{site_id}] downloading + mosaicking …", flush=True)
download_and_mosaic_aster(site, paths, interactive_login=False)

print(f"[{site_id}] regenerating figures …", flush=True)
run_site(site, repo_root, download=False, skip_figures=False)

print(f"[{site_id}] done", flush=True)
"""

_FIGURES_ONLY_WORKER = """
import json
import sys
from pathlib import Path

site_id   = sys.argv[1]
repo_root = Path(sys.argv[2])
# sys.argv[3] is an optional JSON-encoded global_limits dict
global_limits = json.loads(sys.argv[3]) if len(sys.argv) > 3 else None
# JSON gives lists; convert to tuples expected by the function
if global_limits:
    global_limits = {k: tuple(v) for k, v in global_limits.items()}

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.pipeline import run_site

site = load_site_by_id(site_id, repo_root / "sites")
run_site(site, repo_root, download=False, skip_figures=False, global_limits=global_limits)
print(f"[{site_id}] done", flush=True)
"""


def sites_with_data() -> list[str]:
    """Site IDs that have at least one TIR B10 file on disk."""
    found = []
    for site_dir in sorted(DATA_DIR.iterdir()):
        aster_dir = site_dir / "aster"
        if any(aster_dir.glob("*_TIR_B10.tif")):
            found.append(site_dir.name)
    return found


def _launch(sid: str, worker_script: str, extra_args: list[str] | None = None) -> tuple[str, int, str]:
    """Run one site in a subprocess; return (site_id, returncode, output)."""
    cmd = [sys.executable, "-c", worker_script, sid, str(REPO_ROOT)]
    if extra_args:
        cmd.extend(extra_args)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    combined = (proc.stdout + proc.stderr).strip()
    return sid, proc.returncode, combined


def run_parallel(
    site_ids: list[str],
    worker_script: str,
    max_workers: int,
    extra_args: list[str] | None = None,
) -> None:
    t0 = time.monotonic()
    failed: list[str] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_launch, sid, worker_script, extra_args): sid
            for sid in site_ids
        }
        for fut in as_completed(futures):
            sid, rc, output = fut.result()
            elapsed = time.monotonic() - t0
            status = "OK" if rc == 0 else "FAILED"
            print(f"[{elapsed:5.0f}s] {sid}: {status}")
            if rc != 0:
                failed.append(sid)
                for line in output.splitlines()[-20:]:
                    print(f"       {line}")

    print(
        f"\nFinished in {time.monotonic() - t0:.0f}s — "
        f"{len(site_ids) - len(failed)}/{len(site_ids)} succeeded."
    )
    if failed:
        print(f"Failed: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    import json
    from pathlib import Path as _Path

    args = sys.argv[1:]
    remosaic = "--remosaic" in args
    no_global = "--no-global-limits" in args
    site_args = [a for a in args if not a.startswith("--")]

    site_ids = site_args if site_args else sites_with_data()
    if not site_ids:
        print("No sites with ASTER data found.", file=sys.stderr)
        sys.exit(1)

    mode = "re-download + re-mosaic + figures" if remosaic else "figures only (cached data)"
    print(f"Mode     : {mode}")
    print(f"Sites    : {', '.join(site_ids)}")
    max_w = MAX_PARALLEL_DOWNLOADS if remosaic else len(site_ids)
    print(f"Workers  : {max_w}\n")

    worker = _REMOSAIC_WORKER if remosaic else _FIGURES_ONLY_WORKER
    extra_args: list[str] | None = None

    if not remosaic and not no_global:
        # Compute cross-site colorbar limits before spawning workers so every
        # site's Figure 01 uses the same vmin/vmax for comparable colorbars.
        print("Computing cross-site global band ratio limits …", flush=True)
        sys.path.insert(0, str(REPO_ROOT / "src"))
        from critical_minerals_aster.pipeline import compute_global_limits  # noqa: E402

        gl = compute_global_limits(site_ids, REPO_ROOT)
        print(f"  silica    : {gl['silica'][0]:.4f} – {gl['silica'][1]:.4f}")
        print(f"  carbonate : {gl['carbonate'][0]:.4f} – {gl['carbonate'][1]:.4f}")
        print(f"  mafic     : {gl['mafic'][0]:.4f} – {gl['mafic'][1]:.4f}\n")
        extra_args = [json.dumps({k: list(v) for k, v in gl.items()})]

    run_parallel(site_ids, worker, max_workers=max_w, extra_args=extra_args)
