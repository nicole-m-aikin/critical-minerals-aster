#!/usr/bin/env python
"""Goldfield AVIRIS-5 mineral-index pilot.

What this does, in order:
  1. Pick AVIRIS-5 scenes that together paint the ASTER study bbox
     (greedy set-cover, seeded with town/district), same idea as the
     ASTER --mosaic path.
  2. Download (or reuse) those ~2 GB surface-reflectance cubes.
  3. Read only the wavelengths we need (~2120–2270 nm), not all 424 bands.
  4. Measure how deep the spectrum dips at clay/alunite wavelengths.
  5. Split that into three maps:
       Al-OH     — any 2.20 um dip (clay *or* alunite *or* mica)
       kaolinite — *both* 2.16 and 2.21 um dips (the kaolinite pair)
       alunite   — *both* 2.17 and 2.21 um dips (the alunite pair)

Goldfield is a high-sulfidation epithermal district.  The textbook surface
mineral is alunite.  If the alunite panel is blank, the chain is wrong.

Usage:
    conda run -n aster-minerals python scripts/gemx_goldfield_pilot.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.aviris import (
    download_rfl_ort,
    mosaic_to_bbox,
    read_rfl_subset,
    select_granules_painting_bbox,
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

SITE_ID = "goldfield"
CACHE = ROOT / "data" / "gemx" / "cmr_cache" / "goldfield__aviris5_l2a_rfl.json"
OUT_DIR = ROOT / "data" / "gemx" / "aviris" / "goldfield"
FIG_DIR = ROOT / "figures"
RESULTS = ROOT / "results"
# Same WGS84 points as the ASTER Goldfield figures / GNIS.
# Town: Mindat / GNIS populated place. District: GNIS mining district.
TOWN = (-117.23563, 37.70854)
DISTRICT = (-117.14341, 37.71938)

# One shared continuum for the whole Al-OH window.  Shoulders sit between
# the 1.9 um water band and the 2.3 um Mg-OH band.
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


def _depth_at(stack, actual_w, feature_nm):
    i_lo, i_feat, i_hi = 0, None, len(actual_w) - 1
    i_feat = int(np.argmin(np.abs(actual_w - feature_nm)))
    return continuum_removed_band_depth(
        stack[i_lo], stack[i_feat], stack[i_hi],
        float(actual_w[i_lo]), float(actual_w[i_feat]), float(actual_w[i_hi]),
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


def main() -> None:
    site = load_site_by_id(SITE_ID, ROOT / "sites")
    lon = (site.bbox_wgs84[0] + site.bbox_wgs84[2]) / 2
    lat = (site.bbox_wgs84[1] + site.bbox_wgs84[3]) / 2
    town_lon, town_lat = TOWN
    granules = select_granules_painting_bbox(
        CACHE, site.bbox_wgs84, must_contain=(TOWN, DISTRICT), min_gain=0.02,
    )
    if not granules:
        raise FileNotFoundError("no AVIRIS-5 footprint intersects the ASTER study bbox")
    print(f"Granules painting the ASTER bbox ({len(granules)}):")
    for g in granules:
        print(f"  {g}")

    per_scene = {k: [] for k in ("aloh", "kaolinite", "alunite")}
    transforms = []
    nc_paths = []
    wavelengths = None
    features = None
    crs = None
    used = []
    for granule_ur in granules:
        try:
            nc_path = download_rfl_ort(granule_ur, OUT_DIR)
            stack, actual_w, transform, crs = read_rfl_subset(nc_path, TARGETS)
        except Exception as exc:
            print(f"  skip {granule_ur}: {exc}", flush=True)
            continue
        used.append(granule_ur)
        nc_paths.append(str(nc_path.relative_to(ROOT)))
        print(f"{granule_ur} wavelengths (nm):",
              ", ".join(f"{w:.1f}" for w in actual_w))
        panels, features = _indices_for_cube(stack, actual_w)
        for key, arr in panels.items():
            per_scene[key].append(arr)
        transforms.append(transform)
        wavelengths = [float(x) for x in actual_w]
    if not used:
        raise RuntimeError("no AVIRIS-5 cubes could be read")
    granules = used

    mosaics = {}
    mosaic_t = None
    mosaic_crs = None
    for key, layers in per_scene.items():
        mosaics[key], mosaic_t, mosaic_crs = mosaic_to_bbox(
            layers, transforms, crs, site.bbox_wgs84, resolution=10.0,
        )

    stats = {
        "site_id": SITE_ID,
        "granule_ur": granules,
        "nc_path": nc_paths,
        "wavelengths_nm": wavelengths,
        "features_used_nm": features,
        "crs": "EPSG:32611",
        "town_lon": town_lon,
        "town_lat": town_lat,
        "district_lon": DISTRICT[0],
        "district_lat": DISTRICT[1],
        "bbox_centroid": [lon, lat],
        "site_bbox_wgs84": list(site.bbox_wgs84),
        "mosaic_resolution_m": 10.0,
        "mosaic_bbox_frac": float(np.mean(np.isfinite(next(iter(mosaics.values()))))),
        "indices_site_bbox": {k: _summarize(k, v) for k, v in mosaics.items()},
        "note": (
            "AVIRIS-5 mosaic onto site.bbox_wgs84 / EPSG:32611, the same "
            "window as the ASTER figures. Scenes chosen by greedy bbox "
            "set-cover (ASTER mosaic analogue; flight lines cannot pass a "
            "30% single-scene gate), seeded with town and district. "
            "First finite pixel wins. Not Tetracorder."
        ),
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    stats_path = RESULTS / "gemx_goldfield_pilot.json"
    stats_path.write_text(json.dumps(stats, indent=2))
    print("on ASTER bbox:", json.dumps(stats["indices_site_bbox"], indent=2))

    for key, arr in mosaics.items():
        write_geotiff(OUT_DIR / f"goldfield_{key}.tif", arr, mosaic_t, mosaic_crs)

    fig_path = _save_figure(mosaics, mosaic_t, mosaic_crs, site, stats)
    print(f"Wrote {stats_path.relative_to(ROOT)}")
    print(f"Wrote {fig_path.relative_to(ROOT)}")


def _figure_from_tifs(site) -> Path:
    """Redraw from mosaicked GeoTIFFs so hillshade/outline tweaks skip the cubes."""
    import rasterio

    stats = json.loads((RESULTS / "gemx_goldfield_pilot.json").read_text())
    panels = {}
    transform = crs = None
    for key in ("aloh", "kaolinite", "alunite"):
        with rasterio.open(OUT_DIR / f"goldfield_{key}.tif") as src:
            arr = src.read(1).astype(np.float32)
            nodata = src.nodata
            if nodata is not None:
                arr = np.where(arr == nodata, np.nan, arr)
            panels[key] = arr
            transform = src.transform
            crs = src.crs
    return _save_figure(panels, transform, crs, site, stats)


def _town_xy(crs, lon, lat):
    from pyproj import Transformer

    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    return to_utm.transform(lon, lat)


def _bbox_label(bbox) -> str:
    west, south, east, north = bbox
    return (
        f"{south:.2f}°–{north:.2f}°N · "
        f"{abs(west):.2f}°–{abs(east):.2f}°W"
    )


def _load_site_basemap():
    """Same Esri imagery and UTM window as the ASTER Goldfield figures."""
    import rasterio
    from rasterio.plot import plotting_extent

    path = ROOT / "data" / "basemap" / "goldfield" / "rgb.tif"
    with rasterio.open(path) as src:
        rgb = np.transpose(src.read([1, 2, 3]), (1, 2, 0))
        extent_m = plotting_extent(src)
    return rgb, list(extent_m)


def _hillshade_overlay(hillshade: np.ndarray, strength: float = 0.28) -> np.ndarray:
    """Transparent shadow overlay; the satellite/mosaic colors stay underneath."""
    hs = np.asarray(hillshade, dtype=np.float32)
    valid = np.isfinite(hs)
    hs = np.where(valid, np.clip(hs, 0.0, 1.0), 1.0)
    rgba = np.zeros(hs.shape + (4,), dtype=np.float32)
    rgba[..., 3] = np.where(valid, (1.0 - hs) * strength, 0.0)
    return rgba


def _load_hillshade(site, crs):
    """NASADEM hillshade on the site bbox at ASTER-like 90 m."""
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
    extent = plotting_extent(
        np.empty(hs_shape), transform=hs_transform,
    )
    return hillshade, list(extent)


def _save_figure(panels: dict, transform, crs, site, stats: dict) -> Path:
    import matplotlib.ticker
    from matplotlib.colors import PowerNorm
    from matplotlib.lines import Line2D
    from matplotlib_scalebar.scalebar import ScaleBar
    from rasterio.plot import plotting_extent

    from critical_minerals_aster.spectral import compute_valid_data_footprint

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
    dist_x, dist_y = _town_xy(crs, stats["district_lon"], stats["district_lat"])
    rgb, bbox_extent = _load_site_basemap()
    hillshade, hs_extent = _load_hillshade(site, crs)
    mosaic_extent = list(plotting_extent(panels["aloh"], transform))
    from rasterio.crs import CRS as RasterioCRS
    footprint = compute_valid_data_footprint(
        panels["aloh"], transform, RasterioCRS.from_user_input(crs),
    )
    cover_pct = stats.get("mosaic_bbox_frac", 0.0) * 100.0

    last_im = None
    for i, (ax, key) in enumerate(zip(axes, ("aloh", "kaolinite", "alunite"))):
        ax.set_facecolor("#d0d0d0")
        ax.imshow(rgb, extent=bbox_extent, origin="upper", zorder=0,
                  interpolation="bilinear")
        data = np.ma.masked_invalid(panels[key])
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
        ax.annotate(
            "Goldfield",
            (town_x, town_y),
            textcoords="offset points", xytext=(-8, 8),
            fontsize=8, color="#111111", ha="right",
        )
        ax.plot(
            dist_x, dist_y, marker="^", color="black", markersize=7,
            markeredgecolor="white", markeredgewidth=0.4, zorder=5,
        )
        ax.annotate(
            "district",
            (dist_x, dist_y),
            textcoords="offset points", xytext=(8, 6),
            fontsize=7, color="#111111", ha="left",
        )
        ax.set_title(titles[key], fontsize=10)
        ax.set_xlabel("Easting (km, EPSG:32611)")
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

    axes[0].set_ylabel("Northing (km, EPSG:32611)")
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
        "Goldfield–Cuprite, NV  ·  AVIRIS-5 mosaic on the ASTER study bbox\n"
        f"{_bbox_label(site.bbox_wgs84)}   ·   {len(granules)} granules   ·   "
        "not a spectral-library ID",
        fontsize=11,
    )
    out = FIG_DIR / "gemx_goldfield_aloh.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


if __name__ == "__main__":
    if "--from-tif" in sys.argv:
        site = load_site_by_id(SITE_ID, ROOT / "sites")
        path = _figure_from_tifs(site)
        print(f"Wrote {path.relative_to(ROOT)}")
    else:
        main()
