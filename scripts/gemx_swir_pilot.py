"""AVIRIS Al-OH / kaolinite / alunite mosaic onto an ASTER study bbox.

Collection-aware: census selects AVIRIS-5 or Classic per site.  Batch-capable
via ``--primary-fdr`` / ``--all-covered`` / ``--sites``.  Not Tetracorder.

Usage:
    conda run -n aster-minerals python scripts/gemx_swir_pilot.py --site marysvale
    conda run -n aster-minerals python scripts/gemx_swir_pilot.py --site marysvale --from-tif
    conda run -n aster-minerals python scripts/gemx_swir_pilot.py --primary-fdr --from-tif
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.aviris import (
    COLLECTION_GSD_M,
    COLLECTION_SENSORS,
    download_swir_product,
    mosaic_to_bbox,
    read_swir_subset,
    select_granules_painting_bbox,
    wavelength_targets_ok,
    write_geotiff,
)
from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.hyperspectral import (
    ALOH_CONTINUUM_NM,
    ALOH_FEATURE_NM,
    ALUNITE_FEATURES_NM,
    KAOLINITE_FEATURES_NM,
    continuum_removed_band_depth,
    doublet_score,
)
from gemx_site_select import (  # noqa: E402
    AVIRIS5_COLLECTION,
    SITE_PINS,
    add_site_select_args,
    aviris5_coverage_by_site,
    census_collection_by_site,
    resolve_pins,
    select_census_collection,
    sites_from_namespace,
)

# Re-export for gemx_swir_sensitivity / correspondence imports.
__all__ = ["SITE_PINS", "resolve_pins", "main", "run_site", "swir_tifs_exist"]

SHOULDERS = ALOH_CONTINUUM_NM
TARGETS = (
    SHOULDERS[0],
    KAOLINITE_FEATURES_NM[0],
    ALUNITE_FEATURES_NM[0],
    ALOH_FEATURE_NM,
    KAOLINITE_FEATURES_NM[1],
    ALUNITE_FEATURES_NM[1],
    SHOULDERS[1],
)


def _site_id_from_argv(default: str = "marysvale") -> str:
    """Legacy single-site argv parse (kept for sensitivity scripts)."""
    if "--site" in sys.argv:
        i = sys.argv.index("--site")
        if i + 1 < len(sys.argv):
            return sys.argv[i + 1]
    return default


def _paths(site_id: str, collection_key: str = AVIRIS5_COLLECTION) -> dict:
    return {
        "collection_key": collection_key,
        "sensor": COLLECTION_SENSORS.get(collection_key, collection_key),
        "cache": ROOT / "data" / "gemx" / "cmr_cache" / f"{site_id}__{collection_key}.json",
        "out": ROOT / "data" / "gemx" / "aviris" / site_id,
        "fig": ROOT / "figures" / f"gemx_{site_id}_aloh.png",
        "stats": ROOT / "results" / f"gemx_{site_id}_pilot.json",
        "provenance": ROOT / "results" / f"gemx_{site_id}_swir_provenance.json",
        "basemap": ROOT / "data" / "basemap" / site_id / "rgb.tif",
    }


def _census_pick(site_id: str) -> tuple[str, float, str] | None:
    picked = census_collection_by_site(ROOT)
    return picked.get(site_id)


SWIR_KEYS = ("aloh", "kaolinite", "alunite")


def swir_tifs_exist(site_id: str, paths: dict | None = None) -> bool:
    """True if Al-OH / kaolinite / alunite GeoTIFFs are already on disk."""
    p = paths or _paths(site_id)
    return all((p["out"] / f"{site_id}_{key}.tif").is_file() for key in SWIR_KEYS)


def _depth_at(stack, actual_w, feature_nm):
    i_feat = int(np.argmin(np.abs(actual_w - feature_nm)))
    return continuum_removed_band_depth(
        stack[0], stack[i_feat], stack[-1],
        float(actual_w[0]), float(actual_w[i_feat]), float(actual_w[-1]),
    ), float(actual_w[i_feat])


def _summarize(name, arr):
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return {"name": name, "n_valid": 0}
    return {
        "name": name,
        "n_valid": int(finite.size),
        "p50": float(np.nanpercentile(finite, 50)),
        "p90": float(np.nanpercentile(finite, 90)),
        "p99": float(np.nanpercentile(finite, 99)),
        "frac_gt_0_05": float(np.mean(finite > 0.05)),
        "frac_gt_0_10": float(np.mean(finite > 0.10)),
    }


def _indices_for_cube(stack, actual_w):
    aloh, w_aloh = _depth_at(stack, actual_w, ALOH_FEATURE_NM)
    k2165, w_k1 = _depth_at(stack, actual_w, KAOLINITE_FEATURES_NM[0])
    k2206, w_k2 = _depth_at(stack, actual_w, KAOLINITE_FEATURES_NM[1])
    a2175, w_a1 = _depth_at(stack, actual_w, ALUNITE_FEATURES_NM[0])
    a2210, w_a2 = _depth_at(stack, actual_w, ALUNITE_FEATURES_NM[1])
    return {
        "aloh": aloh,
        "kaolinite": doublet_score(k2165, k2206),
        "alunite": doublet_score(a2175, a2210),
    }, {
        "aloh": w_aloh,
        "kaolinite": [w_k1, w_k2],
        "alunite": [w_a1, w_a2],
    }


def _town_xy(crs, lon, lat):
    from pyproj import Transformer

    return Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)


def _bbox_label(bbox) -> str:
    west, south, east, north = bbox
    return (
        f"{south:.2f}°–{north:.2f}°N · "
        f"{abs(west):.2f}°–{abs(east):.2f}°W"
    )


def _load_site_basemap(site_id: str):
    import rasterio
    from rasterio.plot import plotting_extent

    path = ROOT / "data" / "basemap" / site_id / "rgb.tif"
    with rasterio.open(path) as src:
        rgb = np.transpose(src.read([1, 2, 3]), (1, 2, 0))
        extent_m = plotting_extent(src)
    return rgb, list(extent_m)


def _hillshade_overlay(hillshade: np.ndarray, strength: float = 0.28) -> np.ndarray:
    hs = np.asarray(hillshade, dtype=np.float32)
    valid = np.isfinite(hs)
    hs = np.where(valid, np.clip(hs, 0.0, 1.0), 1.0)
    rgba = np.zeros(hs.shape + (4,), dtype=np.float32)
    rgba[..., 3] = np.where(valid, (1.0 - hs) * strength, 0.0)
    return rgba


def _load_hillshade(site, crs):
    from rasterio.crs import CRS
    from rasterio.transform import from_bounds
    from rasterio.warp import transform_bounds
    from rasterio.plot import plotting_extent

    from critical_minerals_aster.paths import site_paths_for
    from critical_minerals_aster.terrain import compute_hillshade_for_site

    dst_crs = CRS.from_user_input(crs)
    west, south, east, north = transform_bounds(
        "EPSG:4326", dst_crs, *site.bbox_wgs84, densify_pts=21,
    )
    px = 90.0
    cols = max(1, int(round((east - west) / px)))
    rows = max(1, int(round((north - south) / px)))
    hs_t = from_bounds(west, south, east, north, cols, rows)
    paths = site_paths_for(site, ROOT)
    result = compute_hillshade_for_site(site, paths, hs_t, (rows, cols), dst_crs)
    if result is None:
        return None, None
    hillshade, hs_transform, hs_shape = result
    extent = plotting_extent(np.empty(hs_shape), transform=hs_transform)
    return hillshade, list(extent)


def _epsg_label(crs) -> str:
    from rasterio.crs import CRS

    obj = CRS.from_user_input(crs)
    if obj.to_epsg():
        return f"EPSG:{obj.to_epsg()}"
    return str(obj)


def _save_figure(panels, transform, crs, site, stats: dict, fig_path: Path) -> Path:
    import matplotlib.ticker
    from matplotlib.colors import PowerNorm
    from matplotlib.lines import Line2D
    from matplotlib_scalebar.scalebar import ScaleBar
    from rasterio.crs import CRS as RasterioCRS
    from rasterio.plot import plotting_extent

    from critical_minerals_aster.spectral import compute_valid_data_footprint

    pins = resolve_pins(site, sites_dir=ROOT / "sites")
    titles = {
        "aloh": "a  Any Al–OH clay\n(one dip at 2.20 μm)",
        "kaolinite": "b  Kaolinite-like\n(both 2.16 and 2.21 μm)",
        "alunite": "c  Alunite-like\n(both 2.17 and 2.21 μm)",
    }
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 6.8), constrained_layout=True)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad(alpha=0.0)
    norm = PowerNorm(gamma=0.55, vmin=0.0, vmax=0.12)
    town_x, town_y = _town_xy(crs, stats["town_lon"], stats["town_lat"])
    tgt_x, tgt_y = _town_xy(crs, stats["target_lon"], stats["target_lat"])
    rgb, bbox_extent = _load_site_basemap(site.id)
    hillshade, hs_extent = _load_hillshade(site, crs)
    mosaic_extent = list(plotting_extent(panels["aloh"], transform))
    footprint = compute_valid_data_footprint(
        panels["aloh"], transform, RasterioCRS.from_user_input(crs),
    )
    cover_pct = stats.get("mosaic_bbox_frac", 0.0) * 100.0
    epsg = _epsg_label(crs)

    last_im = None
    for i, (ax, key) in enumerate(zip(axes, ("aloh", "kaolinite", "alunite"))):
        ax.set_facecolor("#d0d0d0")
        ax.imshow(rgb, extent=bbox_extent, origin="upper", zorder=0,
                  interpolation="bilinear")
        data = np.ma.masked_invalid(panels[key])
        if key == "aloh":
            # Negative depth is a peak, not a clay dip; hide it so the
            # panel is not a cream wash (Marysvale vegetation/continuum).
            data = np.ma.masked_where(data <= 0, data)
        last_im = ax.imshow(
            data, extent=mosaic_extent, cmap=cmap, norm=norm,
            interpolation="nearest", origin="upper", zorder=1,
        )
        if hillshade is not None:
            ax.imshow(
                _hillshade_overlay(hillshade),
                extent=hs_extent, origin="upper", zorder=3,
                interpolation="bilinear",
            )
        if footprint is not None and len(footprint):
            footprint.plot(
                ax=ax, facecolor="none", edgecolor="#1a252f",
                linewidth=2.2, linestyle="--", alpha=0.85, zorder=4,
            )
        ax.set_xlim(bbox_extent[0], bbox_extent[1])
        ax.set_ylim(bbox_extent[2], bbox_extent[3])
        ax.set_aspect("equal")
        ax.plot(
            town_x, town_y, marker="*", color="black", markersize=11,
            markeredgecolor="white", markeredgewidth=0.5, zorder=5,
        )
        town_off = pins.get("town_offset", (-8, 8))
        ax.annotate(
            pins["town_label"], (town_x, town_y),
            textcoords="offset points", xytext=town_off,
            fontsize=8, color="#111111", ha="right",
        )
        ax.plot(
            tgt_x, tgt_y, marker="^", color="black", markersize=7,
            markeredgecolor="white", markeredgewidth=0.4, zorder=5,
        )
        tgt_off = pins.get("target_offset", (8, 6))
        ax.annotate(
            pins["target_label"], (tgt_x, tgt_y),
            textcoords="offset points", xytext=tgt_off,
            fontsize=7, color="#111111", ha="left",
        )
        ax.set_title(titles[key], fontsize=10)
        ax.set_xlabel(f"Easting (km, {epsg})")
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.tick_params(labelsize=7)
        s = stats["indices_site_bbox"][key]
        ax.text(
            0.98, 0.02,
            f"in box  p90 = {s.get('p90', float('nan')):.3f}\n"
            f">0.10: {s.get('frac_gt_0_10', 0.0) * 100:.1f}%",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color="#333333",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=2),
        )
        if i == 0:
            ax.add_artist(ScaleBar(
                1, "m", length_fraction=0.15, location="lower left",
                box_alpha=0.7, font_properties={"size": 8},
            ))
            ax.annotate(
                "N\n↑", xy=(0.04, 0.96), xycoords="axes fraction",
                ha="center", va="top", fontsize=11, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          alpha=0.8, edgecolor="gray"),
            )
            ax.text(
                0.02, 0.14,
                f"SWIR data covers {cover_pct:.0f}% of site bbox",
                transform=ax.transAxes, fontsize=8, va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          alpha=0.85, edgecolor="none"),
                zorder=10,
            )

    axes[0].set_ylabel(f"Northing (km, {epsg})")
    fig.legend(
        handles=[
            Line2D([0], [0], color="#1a252f", linestyle="--", linewidth=2.2,
                   label="SWIR data boundary"),
        ],
        loc="lower center", ncol=1, frameon=False, fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
    )
    cbar = fig.colorbar(
        last_im, ax=axes, fraction=0.02, pad=0.02,
        label="Band depth  (1 − R / continuum)",
    )
    cbar.set_ticks([0.00, 0.02, 0.04, 0.08, 0.12])
    granules = stats["granule_ur"]
    if isinstance(granules, str):
        granules = [granules]
    fig.suptitle(
        f"{site.name}  ·  AVIRIS-5 mosaic on the ASTER study bbox\n"
        f"{_bbox_label(site.bbox_wgs84)}   ·   {len(granules)} granules   ·   "
        "not a spectral-library ID",
        fontsize=11,
    )
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def _figure_from_tifs(site, paths: dict) -> Path:
    import rasterio

    stats = json.loads(paths["stats"].read_text())
    panels = {}
    transform = crs = None
    for key in ("aloh", "kaolinite", "alunite"):
        with rasterio.open(paths["out"] / f"{site.id}_{key}.tif") as src:
            arr = src.read(1).astype(np.float32)
            nodata = src.nodata
            if nodata is not None:
                arr = np.where(arr == nodata, np.nan, arr)
            panels[key] = arr
            transform = src.transform
            crs = src.crs
    return _save_figure(panels, transform, crs, site, stats, paths["fig"])


def run_site(
    site_id: str,
    *,
    from_tif: bool = False,
    skip_if_no_coverage: bool = True,
) -> dict:
    """Run band-depth mosaic (or figure-from-tif) for one site.

    Returns a status dict suitable for batch CSV rows.  Does not raise on
    missing coverage when ``skip_if_no_coverage`` is True.
    """
    site = load_site_by_id(site_id, ROOT / "sites")
    pins = resolve_pins(site, sites_dir=ROOT / "sites")
    pick = _census_pick(site_id)
    if pick is None:
        row = {
            "site_id": site_id,
            "status": "no_coverage",
            "swir_coverage": 0.0,
            "collection": "",
            "sensor": "",
            "pins_source": pins.get("pins_source", "site_pins"),
        }
        if skip_if_no_coverage:
            print(f"[{site_id}] skip: no full-SWIR coverage", flush=True)
            return row
        raise FileNotFoundError(f"{site_id}: no full-SWIR census collection")
    collection_key, cov, sensor = pick
    paths = _paths(site_id, collection_key)

    row = {
        "site_id": site_id,
        "status": "ok",
        "swir_coverage": cov,
        "collection": collection_key,
        "sensor": sensor,
        "pins_source": pins.get("pins_source", "site_pins"),
    }

    have_tifs = swir_tifs_exist(site_id, paths)
    # Existing AVIRIS-5 TIFFs are only resumable when census selected AVIRIS-5.
    if have_tifs and collection_key != AVIRIS5_COLLECTION:
        # Allow resume for Classic only when a provenance file matches.
        prov = paths["provenance"]
        if prov.is_file():
            import json as _json
            meta = _json.loads(prov.read_text())
            if meta.get("collection") == collection_key:
                pass  # OK to resume
            else:
                have_tifs = False
        else:
            # Legacy AVIRIS-5 TIFFs under a Classic census pick need reprocess.
            have_tifs = False
    if skip_if_no_coverage and cov <= 0 and not have_tifs:
        row["status"] = "no_coverage"
        print(f"[{site_id}] skip: no SWIR coverage", flush=True)
        return row

    # Resume: existing GeoTIFFs skip multi-GB cube download.
    if have_tifs:
        try:
            fig_path = _figure_from_tifs(site, paths)
            print(f"[{site_id}] Wrote {fig_path.relative_to(ROOT)} (from-tif)", flush=True)
            row["status"] = "ok_from_tif"
        except Exception as exc:
            row["status"] = "error"
            row["error"] = str(exc)
            print(f"[{site_id}] error: {exc}", flush=True)
        return row

    if from_tif:
        row["status"] = "insufficient"
        print(f"[{site_id}] --from-tif but GeoTIFFs missing", flush=True)
        return row

    if not paths["cache"].is_file() and cov <= 0:
        row["status"] = "no_coverage"
        print(f"[{site_id}] skip: no CMR cache / coverage", flush=True)
        return row

    try:
        main(site_id, collection_key=collection_key)
        row["status"] = "ok"
    except Exception as exc:
        row["status"] = "error"
        row["error"] = str(exc)
        print(f"[{site_id}] error: {exc}", flush=True)
    return row


def main(site_id: str, collection_key: str | None = None) -> None:
    site = load_site_by_id(site_id, ROOT / "sites")
    pins = resolve_pins(site, sites_dir=ROOT / "sites")
    if collection_key is None:
        pick = _census_pick(site_id)
        if pick is None:
            raise FileNotFoundError(f"{site_id}: no full-SWIR census collection")
        collection_key, cov, sensor = pick
    else:
        sensor = COLLECTION_SENSORS.get(collection_key, collection_key)
        cov = None
    paths = _paths(site_id, collection_key)
    resolution = float(COLLECTION_GSD_M.get(collection_key, 10.0))
    town = pins["town"]
    target = pins["target"]
    if not paths["cache"].is_file():
        raise FileNotFoundError(
            f"missing CMR cache {paths['cache']}; "
            "run gemx_coverage_report first or place the cache JSON"
        )
    granules = select_granules_painting_bbox(
        paths["cache"], site.bbox_wgs84, must_contain=(town, target), min_gain=0.02,
    )
    if not granules:
        raise FileNotFoundError(
            f"no {sensor} footprint intersects the ASTER study bbox"
        )
    print(f"{sensor} granules painting the ASTER bbox ({len(granules)}):")
    for g in granules:
        print(f"  {g}")

    per_scene = {k: [] for k in SWIR_KEYS}
    transforms = []
    product_paths = []
    wavelengths = None
    features = None
    crs = None
    used = []
    for granule_ur in granules:
        try:
            product = download_swir_product(
                granule_ur, paths["out"], collection_key=collection_key,
            )
            if collection_key == AVIRIS5_COLLECTION:
                cube_path = product
                hdr_path = None
            else:
                hdr_path, cube_path = product
            stack, actual_w, transform, crs = read_swir_subset(
                cube_path, TARGETS, collection_key=collection_key, hdr_path=hdr_path,
            )
            if not wavelength_targets_ok(actual_w, TARGETS):
                raise ValueError("insufficient_spectral_sampling")
        except Exception as exc:
            print(f"  skip {granule_ur}: {exc}", flush=True)
            continue
        used.append(granule_ur)
        if collection_key == AVIRIS5_COLLECTION:
            product_paths.append(str(Path(cube_path).relative_to(ROOT)))
        else:
            product_paths.append({
                "hdr": str(Path(hdr_path).relative_to(ROOT)),
                "bin": str(Path(cube_path).relative_to(ROOT)),
            })
        print(f"{granule_ur} wavelengths (nm):",
              ", ".join(f"{w:.1f}" for w in actual_w), flush=True)
        panels, features = _indices_for_cube(stack, actual_w)
        for key, arr in panels.items():
            per_scene[key].append(arr)
        transforms.append(transform)
        wavelengths = [float(x) for x in actual_w]
    if not used:
        raise RuntimeError(f"no {sensor} cubes could be read")

    mosaics = {}
    mosaic_t = mosaic_crs = None
    for key, layers in per_scene.items():
        mosaics[key], mosaic_t, mosaic_crs = mosaic_to_bbox(
            layers, transforms, crs, site.bbox_wgs84, resolution=resolution,
        )

    stats = {
        "site_id": site_id,
        "collection": collection_key,
        "sensor": sensor,
        "swir_coverage": cov,
        "granule_ur": used,
        "product_paths": product_paths,
        "wavelengths_nm": wavelengths,
        "features_used_nm": features,
        "crs": _epsg_label(mosaic_crs),
        "town_lon": town[0],
        "town_lat": town[1],
        "target_lon": target[0],
        "target_lat": target[1],
        "target_label": pins["target_label"],
        "pins_source": pins.get("pins_source", "site_pins"),
        "bbox_centroid": list(bbox_centroid_xy(site)),
        "site_bbox_wgs84": list(site.bbox_wgs84),
        "mosaic_resolution_m": resolution,
        "mosaic_bbox_frac": float(np.mean(np.isfinite(next(iter(mosaics.values()))))),
        "indices_site_bbox": {k: _summarize(k, v) for k, v in mosaics.items()},
        "note": (
            f"{sensor} mosaic onto site.bbox_wgs84, same window as the ASTER "
            "figures. Greedy bbox set-cover seeded with town and the science "
            "target. First finite pixel wins. Not Tetracorder."
        ),
    }
    paths["stats"].parent.mkdir(parents=True, exist_ok=True)
    paths["stats"].write_text(json.dumps(stats, indent=2))
    paths["provenance"].write_text(json.dumps({
        "site_id": site_id,
        "collection": collection_key,
        "sensor": sensor,
        "granule_ur": used,
        "mosaic_resolution_m": resolution,
        "wavelengths_nm": wavelengths,
        "features_used_nm": features,
    }, indent=2))
    print("on ASTER bbox:", json.dumps(stats["indices_site_bbox"], indent=2))
    for key, arr in mosaics.items():
        write_geotiff(paths["out"] / f"{site_id}_{key}.tif", arr, mosaic_t, mosaic_crs)
    fig_path = _save_figure(mosaics, mosaic_t, mosaic_crs, site, stats, paths["fig"])
    print(f"Wrote {paths['stats'].relative_to(ROOT)}")
    print(f"Wrote {fig_path.relative_to(ROOT)}")


def bbox_centroid_xy(site) -> tuple[float, float]:
    from gemx_site_select import bbox_centroid

    return bbox_centroid(site.bbox_wgs84)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="AVIRIS-5 Al-OH / kaolinite / alunite mosaic on ASTER bbox",
    )
    add_site_select_args(p)
    p.add_argument(
        "--from-tif",
        action="store_true",
        help="Rebuild figures from existing GeoTIFFs; skip cube download",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    # Backward-compatible default: single marysvale when no selection flags.
    if not any([args.site, args.sites, args.primary_fdr, args.all_covered]):
        site_ids = ["marysvale"]
    else:
        site_ids = sites_from_namespace(args, ROOT)
    for sid in site_ids:
        run_site(sid, from_tif=args.from_tif)
