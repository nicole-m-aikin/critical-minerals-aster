"""Compare ASTER TIR and AVIRIS-5 SWIR footprints on the Goldfield bbox."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import Point, box

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gemx_goldfield_pilot import (  # noqa: E402
    _hillshade_overlay,
    _load_hillshade,
    _load_site_basemap,
    _town_xy,
)

CRS = "EPSG:32611"
TOWN = (-117.23563, 37.70854)
DISTRICT = (-117.14341, 37.71938)


def _km2(geom) -> float:
    return float(geom.area) / 1e6


def _contains(poly, lon, lat) -> bool:
    return bool(poly.contains(Point(lon, lat)))


def _valid_mask(array: np.ndarray, nodata) -> np.ndarray:
    mask = np.isfinite(array)
    if nodata is not None:
        mask &= array != nodata
    return mask


def _mask_to_geom(mask: np.ndarray, transform, crs):
    from rasterio.features import shapes
    from shapely.geometry import shape
    from shapely.ops import unary_union

    if not np.any(mask):
        from shapely.geometry import GeometryCollection
        return GeometryCollection()
    geoms = [
        shape(geom)
        for geom, val in shapes(mask.astype(np.uint8), transform=transform)
        if val == 1
    ]
    if not geoms:
        from shapely.geometry import GeometryCollection
        return GeometryCollection()
    return unary_union(geoms).simplify(180.0)


def _valid_on_grid(path: Path, valid_fn, dst_t, shape, dst_crs):
    import rasterio
    from rasterio.warp import reproject, Resampling

    rows, cols = shape
    with rasterio.open(path) as src:
        src_valid = valid_fn(src.read(1).astype(np.float32), src.nodata).astype(np.uint8)
        dst = np.zeros((rows, cols), dtype=np.uint8)
        reproject(
            src_valid, dst,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=dst_t, dst_crs=dst_crs,
            resampling=Resampling.nearest,
        )
    return dst.astype(bool)


def compute_overlap(site) -> tuple[dict, dict]:
    from rasterio.crs import CRS as RasterioCRS
    from rasterio.transform import from_bounds
    from rasterio.warp import transform_bounds

    bbox_wgs = box(*site.bbox_wgs84)
    prov = json.loads((ROOT / "results" / "goldfield_provenance.json").read_text())

    tir_path = ROOT / "data" / "sites" / "goldfield" / "aster" / "goldfield_mosaic_TIR_B10.tif"
    swir_path = ROOT / "data" / "gemx" / "aviris" / "goldfield" / "goldfield_aloh.tif"

    def tir_valid(a, nodata):
        mask = _valid_mask(a, nodata) & (a > 0)
        return mask

    def swir_valid(a, nodata):
        return _valid_mask(a, nodata)

    bbox_u = gpd.GeoSeries([bbox_wgs], crs="EPSG:4326").to_crs(CRS).iloc[0]

    west, south, east, north = transform_bounds(
        "EPSG:4326", CRS, *site.bbox_wgs84, densify_pts=21,
    )
    px = 90.0
    cols = max(1, int(round((east - west) / px)))
    rows = max(1, int(round((north - south) / px)))
    dst_t = from_bounds(west, south, east, north, cols, rows)
    dst_crs = RasterioCRS.from_epsg(32611)
    tir_r = _valid_on_grid(tir_path, tir_valid, dst_t, (rows, cols), dst_crs)
    swir_r = _valid_on_grid(swir_path, swir_valid, dst_t, (rows, cols), dst_crs)
    n = float(tir_r.size)
    pix = {
        "tir_pct_bbox": 100.0 * tir_r.mean(),
        "swir_pct_bbox": 100.0 * swir_r.mean(),
        "both_pct_bbox": 100.0 * np.mean(tir_r & swir_r),
        "tir_only_pct_bbox": 100.0 * np.mean(tir_r & ~swir_r),
        "swir_only_pct_bbox": 100.0 * np.mean(swir_r & ~tir_r),
        "neither_pct_bbox": 100.0 * np.mean(~tir_r & ~swir_r),
        "union_pct_bbox": 100.0 * np.mean(tir_r | swir_r),
        "swir_of_tir_pct": 100.0 * np.mean(tir_r & swir_r) / max(tir_r.mean(), 1e-9),
        "tir_of_swir_pct": 100.0 * np.mean(tir_r & swir_r) / max(swir_r.mean(), 1e-9),
        "n_pixels_90m": int(n),
    }

    # Draw from the same 90 m masks as the percentages so holes stay holes.
    tir_g = _mask_to_geom(tir_r, dst_t, dst_crs)
    swir_g = _mask_to_geom(swir_r, dst_t, dst_crs)
    both = _mask_to_geom(tir_r & swir_r, dst_t, dst_crs)
    tir_only = _mask_to_geom(tir_r & ~swir_r, dst_t, dst_crs)
    swir_only = _mask_to_geom(swir_r & ~tir_r, dst_t, dst_crs)
    neither = _mask_to_geom(~tir_r & ~swir_r, dst_t, dst_crs)
    union = _mask_to_geom(tir_r | swir_r, dst_t, dst_crs)

    bbox_km2 = _km2(bbox_u)
    layers = {
        "bbox": bbox_u,
        "tir": tir_g,
        "swir": swir_g,
        "both": both,
        "tir_only": tir_only,
        "swir_only": swir_only,
        "neither": neither,
        "union": union,
    }

    tir_wgs = gpd.GeoSeries([tir_g], crs=CRS).to_crs("EPSG:4326").iloc[0]
    swir_wgs = gpd.GeoSeries([swir_g], crs=CRS).to_crs("EPSG:4326").iloc[0]
    town_in = {"tir": _contains(tir_wgs, *TOWN), "swir": _contains(swir_wgs, *TOWN)}
    dist_in = {"tir": _contains(tir_wgs, *DISTRICT), "swir": _contains(swir_wgs, *DISTRICT)}

    stats = {
        "site_id": site.id,
        "method": "valid pixels warped to 90 m site-bbox grid, EPSG:32611",
        "bbox_km2": round(bbox_km2, 2),
        "tir_km2": round(_km2(tir_g), 2),
        "swir_km2": round(_km2(swir_g), 2),
        "both_km2": round(_km2(both), 2),
        "tir_only_km2": round(_km2(tir_only), 2),
        "swir_only_km2": round(_km2(swir_only), 2),
        "neither_km2": round(_km2(neither), 2),
        "union_km2": round(_km2(union), 2),
        **{k: round(v, 1) if isinstance(v, float) else v for k, v in pix.items()},
        "aster_granule_id": prov.get("granule_id"),
        "swir_n_granules": 10,
        "town_in_tir": town_in["tir"],
        "town_in_swir": town_in["swir"],
        "district_in_tir": dist_in["tir"],
        "district_in_swir": dist_in["swir"],
    }
    return stats, layers


def _gdf(geom, crs=CRS) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(geometry=[geom], crs=crs)


def save_figure(site, stats: dict, layers: dict) -> Path:
    rgb, bbox_extent = _load_site_basemap()
    hillshade, hs_extent = _load_hillshade(site, CRS)
    town_x, town_y = _town_xy(CRS, *TOWN)
    dist_x, dist_y = _town_xy(CRS, *DISTRICT)

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 6.6), constrained_layout=True)
    titles = (
        f"a  ASTER TIR\n{stats['tir_pct_bbox']:.0f}% of bbox · 1 granule",
        f"b  AVIRIS-5 SWIR\n{stats['swir_pct_bbox']:.0f}% of bbox · 10 granules",
        f"c  Overlap\n{stats['both_pct_bbox']:.0f}% both · {stats['union_pct_bbox']:.0f}% either",
    )
    fills = (
        [("#c0392b", layers["tir"], 0.28)],
        [("#1a7a6d", layers["swir"], 0.28)],
        [
            ("#6c3483", layers["both"], 0.32),
            ("#c0392b", layers["tir_only"], 0.28),
            ("#1a7a6d", layers["swir_only"], 0.28),
        ],
    )
    outlines = (
        [("#1a252f", layers["tir"])],
        [("#1a252f", layers["swir"])],
        [("#1a252f", layers["tir"]), ("#0b5345", layers["swir"])],
    )

    for i, ax in enumerate(axes):
        ax.set_facecolor("#d0d0d0")
        ax.imshow(rgb, extent=bbox_extent, origin="upper", zorder=0,
                  interpolation="bilinear")
        if hillshade is not None:
            ax.imshow(
                _hillshade_overlay(hillshade),
                extent=hs_extent, origin="upper", zorder=1,
                interpolation="bilinear",
            )
        for color, geom, alpha in fills[i]:
            if geom.is_empty:
                continue
            _gdf(geom).plot(
                ax=ax, facecolor=color, edgecolor="none", alpha=alpha, zorder=2,
            )
        for color, geom in outlines[i]:
            if geom.is_empty:
                continue
            _gdf(geom).plot(
                ax=ax, facecolor="none", edgecolor=color,
                linewidth=2.0, linestyle="--", alpha=0.9, zorder=4,
            )
        ax.plot(
            town_x, town_y, marker="*", color="black", markersize=11,
            markeredgecolor="white", markeredgewidth=0.5, zorder=5,
        )
        ax.annotate(
            "Goldfield", (town_x, town_y),
            textcoords="offset points", xytext=(-8, 8),
            fontsize=8, color="#111111", ha="right",
        )
        ax.plot(
            dist_x, dist_y, marker="^", color="black", markersize=7,
            markeredgecolor="white", markeredgewidth=0.4, zorder=5,
        )
        ax.set_xlim(bbox_extent[0], bbox_extent[1])
        ax.set_ylim(bbox_extent[2], bbox_extent[3])
        ax.set_aspect("equal")
        ax.set_title(titles[i], fontsize=10)
        ax.set_xlabel("Easting (km, EPSG:32611)")
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.set_ylabel("Northing (km, EPSG:32611)")
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

    fig.legend(
        handles=[
            Patch(facecolor="#c0392b", edgecolor="none", alpha=0.45, label="TIR only"),
            Patch(facecolor="#1a7a6d", edgecolor="none", alpha=0.45, label="SWIR only"),
            Patch(facecolor="#6c3483", edgecolor="none", alpha=0.55, label="Both"),
            Line2D([0], [0], color="#1a252f", linestyle="--", linewidth=2.0,
                   label="TIR boundary"),
            Line2D([0], [0], color="#0b5345", linestyle="--", linewidth=2.0,
                   label="SWIR boundary"),
        ],
        loc="lower center", ncol=5, frameon=False, fontsize=8,
        bbox_to_anchor=(0.5, -0.04),
    )
    fig.suptitle(
        "Goldfield–Cuprite coverage · same ASTER study bbox · "
        "valid-pixel footprints, not mineral ID",
        fontsize=11,
    )
    out = ROOT / "figures" / "gemx_goldfield_coverage_compare.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    from critical_minerals_aster.config import load_site_by_id

    site = load_site_by_id("goldfield", ROOT / "sites")
    stats, layers = compute_overlap(site)
    stats_path = ROOT / "results" / "goldfield_coverage_compare.json"
    # geometries are not JSON-serializable; keep numbers only
    stats_path.write_text(json.dumps(stats, indent=2))
    fig_path = save_figure(site, stats, layers)
    print(json.dumps(stats, indent=2))
    print(f"Wrote {stats_path.relative_to(ROOT)}")
    print(f"Wrote {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
