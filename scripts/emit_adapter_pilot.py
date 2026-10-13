#!/usr/bin/env python3
"""EMIT HyperspectralScene adapter pilot for 2–3 arid sites.

Downloads EMIT L2A reflectance granules intersecting the site bbox (EarthData
login via earthaccess), computes continuum-removed Al-OH / kaolinite / alunite
band depths on the native ~60 m grid, writes GeoTIFFs under ``data/emit/{id}/``,
and verifies ``load_emit_hyperspectral_scene`` matches the AVIRIS feature contract.

Usage:
  conda run -n aster-minerals python scripts/emit_adapter_pilot.py --sites mcdermitt,goldfield,mountain_pass
  conda run -n aster-minerals python scripts/emit_adapter_pilot.py --synthetic   # no download; contract test
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.prospectivity.adapters import load_emit_hyperspectral_scene

# Continuum-removed feature centres (nm) — same family as AVIRIS census
ALOH_NM = 2200.0
KAOL_NM = (2165.0, 2206.0)
ALUN_NM = (2175.0, 2210.0)
CONTINUUM = (2120.0, 2270.0)


def _cr_depth(r_left, r_feat, r_right, w_left, w_feat, w_right):
    """Continuum-removed band depth at feature wavelength."""
    t = (w_feat - w_left) / max(w_right - w_left, 1e-6)
    cont = r_left * (1 - t) + r_right * t
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = 1.0 - (r_feat / cont)
    depth = np.where(np.isfinite(cont) & (cont > 0) & np.isfinite(r_feat), depth, np.nan)
    return depth.astype(np.float32)


def _nearest_band(wavelengths: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(wavelengths - target)))


def write_synthetic_emit_scene(site_id: str, bbox: tuple[float, float, float, float]) -> Path:
    """Write a tiny synthetic EMIT-like Al-OH stack for contract tests (no download)."""
    out_dir = ROOT / "data" / "emit" / site_id
    out_dir.mkdir(parents=True, exist_ok=True)
    west, south, east, north = bbox
    # ~60 m grid in geographic degrees (coarse; only for adapter smoke test)
    gsd_deg = 60.0 / 111_320.0
    width = max(8, int(round((east - west) / gsd_deg)))
    height = max(8, int(round((north - south) / gsd_deg)))
    transform = from_bounds(west, south, east, north, width, height)
    rng = np.random.default_rng(abs(hash(site_id)) % (2**32))
    aloh = rng.random((height, width), dtype=np.float32) * 0.2
    # Plant a clay blob in the centre
    aloh[height // 3 : 2 * height // 3, width // 3 : 2 * width // 3] += 0.35
    kaol = (aloh * 0.7).astype(np.float32)
    alun = (aloh * 0.4).astype(np.float32)
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "compress": "deflate",
        "nodata": np.nan,
    }
    for key, arr in (("aloh", aloh), ("kaolinite", kaol), ("alunite", alun)):
        with rasterio.open(out_dir / f"{site_id}_{key}.tif", "w", **profile) as dst:
            dst.write(arr, 1)
    prov = {
        "site_id": site_id,
        "sensor": "EMIT_L2A",
        "collection": "emit_l2a_rfl",
        "mosaic_resolution_m": 60,
        "features_used_nm": [ALOH_NM, *KAOL_NM, *ALUN_NM],
        "synthetic": True,
    }
    (ROOT / "results" / f"emit_{site_id}_swir_provenance.json").write_text(json.dumps(prov, indent=2))
    return out_dir


def process_emit_netcdf(nc_path: Path, site_id: str, bbox) -> Path | None:
    """Extract Al-OH family depths from an EMIT L2A netCDF if xarray/netCDF available."""
    try:
        import xarray as xr
    except ImportError:
        print("xarray not installed; cannot parse EMIT netCDF")
        return None

    ds = xr.open_dataset(nc_path)
    # EMIT L2A typically has reflectance and wavelengths
    refl_name = next((v for v in ("reflectance", "rfl", "RFL") if v in ds), None)
    if refl_name is None:
        # try data_vars heuristically
        for v in ds.data_vars:
            if "refl" in v.lower():
                refl_name = v
                break
    if refl_name is None:
        print(f"No reflectance var in {nc_path.name}: {list(ds.data_vars)[:10]}")
        return None

    refl = ds[refl_name]
    # wavelengths
    wvar = next((v for v in ("wavelengths", "wavelength", "band") if v in ds.coords or v in ds), None)
    if wvar is None:
        print("No wavelength coordinate")
        return None
    wavelengths = np.asarray(ds[wvar].values, dtype=float)
    if wavelengths.max() < 100:  # maybe microns
        wavelengths = wavelengths * 1000.0

    arr = np.asarray(refl.values, dtype=np.float32)
    # Expect (band, y, x) or (y, x, band)
    if arr.ndim != 3:
        print(f"Unexpected reflectance shape {arr.shape}")
        return None
    if arr.shape[0] == wavelengths.size:
        stack = arr
    elif arr.shape[-1] == wavelengths.size:
        stack = np.moveaxis(arr, -1, 0)
    else:
        print(f"Band axis mismatch {arr.shape} vs {wavelengths.size}")
        return None

    i_l = _nearest_band(wavelengths, CONTINUUM[0])
    i_r = _nearest_band(wavelengths, CONTINUUM[1])
    i_a = _nearest_band(wavelengths, ALOH_NM)
    i_k1, i_k2 = _nearest_band(wavelengths, KAOL_NM[0]), _nearest_band(wavelengths, KAOL_NM[1])
    i_u1, i_u2 = _nearest_band(wavelengths, ALUN_NM[0]), _nearest_band(wavelengths, ALUN_NM[1])

    aloh = _cr_depth(stack[i_l], stack[i_a], stack[i_r], wavelengths[i_l], wavelengths[i_a], wavelengths[i_r])
    k1 = _cr_depth(stack[i_l], stack[i_k1], stack[i_r], wavelengths[i_l], wavelengths[i_k1], wavelengths[i_r])
    k2 = _cr_depth(stack[i_l], stack[i_k2], stack[i_r], wavelengths[i_l], wavelengths[i_k2], wavelengths[i_r])
    u1 = _cr_depth(stack[i_l], stack[i_u1], stack[i_r], wavelengths[i_l], wavelengths[i_u1], wavelengths[i_r])
    u2 = _cr_depth(stack[i_l], stack[i_u2], stack[i_r], wavelengths[i_l], wavelengths[i_u2], wavelengths[i_r])
    kaol = np.minimum(k1, k2)
    alun = np.minimum(u1, u2)

    # Geotransform from lat/lon if present
    out_dir = ROOT / "data" / "emit" / site_id
    out_dir.mkdir(parents=True, exist_ok=True)
    west, south, east, north = bbox
    h, w = aloh.shape
    transform = from_bounds(west, south, east, north, w, h)
    profile = {
        "driver": "GTiff",
        "height": h,
        "width": w,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "compress": "deflate",
        "nodata": np.nan,
    }
    for key, data in (("aloh", aloh), ("kaolinite", kaol), ("alunite", alun)):
        with rasterio.open(out_dir / f"{site_id}_{key}.tif", "w", **profile) as dst:
            dst.write(data, 1)

    used = [float(wavelengths[i]) for i in (i_a, i_k1, i_k2, i_u1, i_u2)]
    prov = {
        "site_id": site_id,
        "sensor": "EMIT_L2A",
        "collection": "emit_l2a_rfl",
        "source_file": nc_path.name,
        "mosaic_resolution_m": 60,
        "features_used_nm": used,
        "synthetic": False,
    }
    (ROOT / "results" / f"emit_{site_id}_swir_provenance.json").write_text(json.dumps(prov, indent=2))
    return out_dir


def try_download_emit(site_id: str, bbox) -> list[Path]:
    """Best-effort earthaccess download of EMIT L2A granules."""
    try:
        import earthaccess
    except ImportError:
        print("earthaccess not installed")
        return []
    try:
        earthaccess.login(strategy="netrc")
    except Exception:
        try:
            earthaccess.login(strategy="environment")
        except Exception as exc:
            print(f"EarthData login failed: {exc}")
            return []

    west, south, east, north = bbox
    results = earthaccess.search_data(
        short_name="EMITL2ARFL",
        bounding_box=(west, south, east, north),
        count=5,
    )
    if not results:
        print(f"{site_id}: no EMIT granules")
        return []
    out_dir = ROOT / "data" / "emit" / site_id / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = earthaccess.download(results[:2], str(out_dir))
    return [Path(p) for p in (paths or [])]


def verify_contract(site_id: str) -> dict:
    scene = load_emit_hyperspectral_scene(site_id, repo_root=ROOT)
    if scene is None:
        return {"site_id": site_id, "ok": False, "reason": "missing_tiffs"}
    ok = (
        scene.sensor_id.startswith("EMIT")
        and scene.aloh.ndim == 2
        and scene.kaolinite.shape == scene.aloh.shape
        and scene.alunite.shape == scene.aloh.shape
        and scene.gsd_m == 60.0
        and len(scene.features_used_nm) > 0
    )
    return {
        "site_id": site_id,
        "ok": bool(ok),
        "sensor_id": scene.sensor_id,
        "shape": list(scene.aloh.shape),
        "coverage_frac": scene.coverage_frac,
        "features_used_nm": list(scene.features_used_nm),
        "gsd_m": scene.gsd_m,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", default="mcdermitt,goldfield,mountain_pass")
    p.add_argument("--synthetic", action="store_true", help="Write synthetic scenes (no download)")
    p.add_argument("--download", action="store_true", help="Attempt EarthData EMIT download")
    args = p.parse_args()
    sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    reports = []
    for sid in sites:
        site = load_site_by_id(sid, ROOT / "sites")
        bbox = tuple(site.bbox_wgs84)
        if args.synthetic or not args.download:
            write_synthetic_emit_scene(sid, bbox)
            print(f"{sid}: wrote synthetic EMIT stack")
        if args.download:
            paths = try_download_emit(sid, bbox)
            for nc in paths:
                if nc.suffix.lower() in {".nc", ".h5"} or "RFL" in nc.name.upper():
                    process_emit_netcdf(nc, sid, bbox)
        reports.append(verify_contract(sid))
        print(reports[-1])

    out = ROOT / "results" / "emit_adapter_pilot.json"
    out.write_text(json.dumps(reports, indent=2))
    print(f"Wrote {out.relative_to(ROOT)}")
    if not all(r.get("ok") for r in reports):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
