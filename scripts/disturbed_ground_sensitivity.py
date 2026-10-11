#!/usr/bin/env python
"""
Verification A — are the strong-anomaly hit rates explained by mine waste?

ASTER TIR strong-anomaly zones at historic districts could be mapping the
spectral contrast of tailings, dumps and pits rather than natural alteration.
If so, deposits "hit" the zones simply because the workings sit on their own
waste.  This script builds a disturbed-ground mask per site from USGS USMIN
mine features (tailings, mine dumps, open pits, disturbed surface polygons,
optionally plus buffered shaft / adit / prospect points; OSM mine-related
polygons as fallback or supplement) and re-runs the critical-mineral hit test.

Variants per site (``variant`` column)
--------------------------------------
    baseline                      -- no mask; reproduces n_crit / hits_crit / p0
    removed_<src>_r<R>            -- mask removed from BOTH zones and TIR footprint
                                     (so p0 stays a consistent area ratio) and
                                     deposits on the mask drop out of the support
    excluded_<src>_r<R>           -- zones / footprint / p0 unchanged; only the
                                     deposits sitting on the mask are dropped
    removed_<src>_halo<H>         -- waste POLYGONS dilated by H m (no point
                                     features), removed from zones + footprint.
                                     Topo-map outlines (1886-2006) can under-map
                                     present-day disturbance (e.g. Holden's
                                     2013-2017 Superfund remediation earthworks).
    Whole connected zone patches touching waste are deliberately NOT removed:
    at several sites one patch holds most of the zone area (Holden: ~80%), so
    that test measures patch connectivity, not waste contamination.
    <src> is "usmin" (USMIN polygons, or OSM polygons where USMIN has none)
    or "usmin+osm" (union with OSM; only when Overpass answered).
    <R> is the point-feature buffer radius in metres (0 = polygons only).

Outputs
-------
    results/disturbed_ground_sensitivity.csv   -- one row per site x variant
    figures/verification_disturbed_ground.png  -- enrichment before/after masking
                                                  + Holden zone / mask / deposit map
    stdout                                     -- summary table, Holden verdict,
                                                  robust sites losing significance

Usage
-----
    conda run -n aster-minerals python scripts/disturbed_ground_sensitivity.py
    conda run -n aster-minerals python scripts/disturbed_ground_sensitivity.py --buffers 0 100 250 --halos 500
    conda run -n aster-minerals python scripts/disturbed_ground_sensitivity.py --no-osm --skip-figure
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from critical_minerals_aster.significance import run_binomial, significant_site_sets
from critical_minerals_aster.usmin import (
    SiteInputs,
    concat_gdfs,
    disturbed_mask,
    feature_type_counts,
    load_site_inputs,
    overlap_fraction,
    points_in_mask,
    subtract_mask,
)

FDR_SITES, ROBUST_SITES = significant_site_sets(REPO_ROOT)

CACHE_DIR = REPO_ROOT / "data" / "usmin"
OUT_CSV = REPO_ROOT / "results" / "disturbed_ground_sensitivity.csv"
OUT_FIG = REPO_ROOT / "figures" / "verification_disturbed_ground.png"
PHASE4_CSV = REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv"


def hit_test(n: int, hits: int, p0: float) -> dict:
    p, expected = run_binomial(hits, n, p0)
    rate = hits / n if n else 0.0
    return {
        "n_crit": n,
        "hits_crit": hits,
        "hit_rate_crit_pct": round(rate * 100, 2),
        "p0_null": round(p0, 4),
        "expected_hits_null": round(expected, 1),
        "enrichment": round(rate / p0, 2) if p0 > 0 else np.nan,
        "p_binomial": p,
        "sig_05": p < 0.05,
    }


def mask_layers(si: SiteInputs, source: str) -> tuple[gpd.GeoSeries, gpd.GeoSeries, str]:
    """(polygons, points, label) in the site CRS for a mask source."""
    usmin_polys = si.usmin_polygons if si.usmin_polygons is not None else None
    have_usmin_polys = usmin_polys is not None and not usmin_polys.empty
    if source == "usmin":
        if have_usmin_polys:
            polys, label = usmin_polys, "usmin"
        elif si.osm_polygons is not None and not si.osm_polygons.empty:
            polys, label = si.osm_polygons, "osm_fallback"
        else:
            polys, label = None, "usmin_points_only"
        pts = si.usmin_points
    else:
        polys = concat_gdfs([usmin_polys, si.osm_polygons], crs=si.crs)
        pts = concat_gdfs([si.usmin_points, si.osm_points], crs=si.crs)
        label = "usmin+osm"
    poly_geoms = polys.to_crs(si.crs).geometry if polys is not None and len(polys) else gpd.GeoSeries([], crs=si.crs)
    pt_geoms = pts.to_crs(si.crs).geometry if pts is not None and len(pts) else gpd.GeoSeries([], crs=si.crs)
    return poly_geoms, pt_geoms, label


def analyse_site(si: SiteInputs, buffers: list[float], halos: list[float]) -> list[dict]:
    dep = si.deposits
    n0, h0 = len(dep), int(dep["inside_zone"].sum())
    fp_km2 = si.footprint.area / 1e6
    zone_km2 = si.zones_union.area / 1e6
    common = {
        "site_id": si.site_id,
        "site_name": si.site_name,
        "robust": si.site_id in ROBUST_SITES,
        "usmin_status": si.usmin_status,
        "osm_status": si.osm_status,
        "n_usmin_points": 0 if si.usmin_points is None else len(si.usmin_points),
        "n_usmin_polygons": 0 if si.usmin_polygons is None else len(si.usmin_polygons),
        "n_osm_points": np.nan if si.osm_points is None else len(si.osm_points),
        "n_osm_polygons": np.nan if si.osm_polygons is None else len(si.osm_polygons),
        "usmin_polygon_types": feature_type_counts(si.usmin_polygons),
    }
    rows = [{
        **common, "variant": "baseline", "mode": "baseline", "mask_source": "",
        "buffer_m": np.nan, "zone_area_km2": round(zone_km2, 3),
        "footprint_area_km2": round(fp_km2, 2), "mask_area_km2": 0.0,
        "pct_zone_on_mask": 0.0, "pct_footprint_on_mask": 0.0, "zone_mask_affinity": np.nan,
        "n_deposits_on_mask": 0, "hits_on_mask": 0,
        **hit_test(n0, h0, si.p0),
    }]

    sources = ["usmin"]
    if si.osm_status == "ok":
        sources.append("usmin+osm")
    for source in sources:
        polys, pts, label = mask_layers(si, source)
        masks = [(f"r{int(r)}", r, 0.0, disturbed_mask(polys, pts, r), ("removed", "excluded"))
                 for r in buffers]
        poly_mask = disturbed_mask(polys, None, 0)
        masks += [(f"halo{int(h)}", 0.0, h, poly_mask.buffer(h) if not poly_mask.is_empty else poly_mask,
                   ("removed",)) for h in halos]
        for suffix, r, h, mask, modes in masks:
            mask = mask.intersection(si.footprint)
            on_mask = points_in_mask(dep.geometry, mask)
            pct_zone = overlap_fraction(si.zones_union, mask) * 100
            pct_fp = overlap_fraction(si.footprint, mask) * 100
            geo = {
                "mask_source": label,
                "buffer_m": r,
                "polygon_halo_m": h,
                "mask_area_km2": round(mask.area / 1e6, 4),
                "pct_zone_on_mask": round(pct_zone, 3),
                "pct_footprint_on_mask": round(pct_fp, 3),
                "zone_mask_affinity": round(pct_zone / pct_fp, 2) if pct_fp > 0 else np.nan,
                "n_deposits_on_mask": int(on_mask.sum()),
                "hits_on_mask": int((on_mask & dep["inside_zone"].values).sum()),
            }
            keep = ~on_mask
            n1, h1 = int(keep.sum()), int((keep & dep["inside_zone"].values).sum())
            tag = f"{source}_{suffix}"
            if "removed" in modes:
                zones_m = subtract_mask(si.zones_union, mask)
                fp_m = subtract_mask(si.footprint, mask)
                p0_m = zones_m.area / fp_m.area if fp_m.area else 0.0
                rows.append({
                    **common, **geo, "variant": f"removed_{tag}", "mode": "removed",
                    "zone_area_km2": round(zones_m.area / 1e6, 3),
                    "footprint_area_km2": round(fp_m.area / 1e6, 2),
                    **hit_test(n1, h1, p0_m),
                })
            if "excluded" in modes:
                rows.append({
                    **common, **geo, "variant": f"excluded_{tag}", "mode": "excluded",
                    "zone_area_km2": round(zone_km2, 3),
                    "footprint_area_km2": round(fp_km2, 2),
                    **hit_test(n1, h1, si.p0),
                })
    return rows


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def make_figure(df: pd.DataFrame, holden: SiteInputs | None, buffers: list[float],
                halos: list[float], out: Path) -> None:
    r_hi = max(buffers)
    h_lo = min(halos) if halos else None
    base = df[df["variant"] == "baseline"].set_index("site_id")
    rem0 = df[df["variant"] == "removed_usmin_r0"].set_index("site_id")
    remh = df[df["variant"] == f"removed_usmin_r{int(r_hi)}"].set_index("site_id")
    order = base.sort_values("enrichment").index.tolist()
    variants = [(rem0, "s", "tab:orange", "USMIN waste polygons removed")]
    if h_lo is not None:
        halo = df[df["variant"] == f"removed_usmin_halo{int(h_lo)}"].set_index("site_id").reindex(base.index)
        variants.append((halo, "^", "tab:green", f"Waste polygons + {int(h_lo)} m halo removed"))
    variants.append((remh, "D", "tab:blue", f"+ {int(r_hi)} m buffers on shafts/adits/pits"))

    fig = plt.figure(figsize=(15, 7.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.22)
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(order))
    for i, s in enumerate(order):
        vals = [base.loc[s, "enrichment"]] + [v.loc[s, "enrichment"] for v, *_ in variants]
        vals = [v for v in vals if pd.notna(v)]
        ax.plot([min(vals), max(vals)], [i, i], color="0.82", lw=1.2, zorder=1)
    ax.scatter(base.loc[order, "enrichment"], y, s=46, color="0.15", label="Baseline (no mask)", zorder=3)
    for v, marker, color, label in variants:
        ax.scatter(v.loc[order, "enrichment"], y, s=48, marker=marker, facecolor="none",
                   edgecolor=color, lw=1.5, label=label, zorder=4)
    for v, _, color, _ in variants:
        for i, s in enumerate(order):
            if pd.notna(v.loc[s, "p_binomial"]) and v.loc[s, "p_binomial"] >= 0.05:
                ax.annotate("p≥0.05", (v.loc[s, "enrichment"], i), xytext=(6, -3),
                            textcoords="offset points", fontsize=7, color=color)
    ax.axvline(1.0, color="0.4", ls="--", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s}{' *' if s in ROBUST_SITES else ''}" for s in order], fontsize=8.5)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("Critical-mineral enrichment (hit rate / site null p0)\n* = fully robust site")
    ax.set_title("A. Enrichment before vs after removing\nmine-disturbed ground (zones + footprint)", fontsize=10.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2, fontsize=7.5, frameon=False)
    ax.set_xlim(left=0)
    ax.grid(axis="x", color="0.92")

    if holden is not None:
        polys, pts, _ = mask_layers(holden, "usmin")
        mask_r0 = disturbed_mask(polys, pts, 0).intersection(holden.footprint)
        mask_rh = disturbed_mask(polys, pts, r_hi).intersection(holden.footprint)
        dep = holden.deposits
        has_osm = holden.osm_polygons is not None and not holden.osm_polygons.empty
        tail = mask_r0.centroid if not mask_r0.is_empty else holden.footprint.centroid
        half = 3500.0
        windows = [
            ("B. Holden, WA — full TIR footprint", holden.footprint.bounds, 12),
            ("C. Holden mine / tailings (7 km window)",
             (tail.x - half, tail.y - half, tail.x + half, tail.y + half), 34),
        ]
        for k, (title, (x0, y0, x1, y1), msize) in enumerate(windows):
            axm = fig.add_subplot(gs[k + 1])
            gpd.GeoSeries([holden.footprint], crs=holden.crs).boundary.plot(ax=axm, color="0.3", lw=0.8)
            gpd.GeoSeries([holden.zones_union], crs=holden.crs).plot(
                ax=axm, color="#f4a259", alpha=0.75, edgecolor="#c0561b", lw=0.3)
            if not mask_rh.is_empty:
                gpd.GeoSeries([mask_rh], crs=holden.crs).plot(ax=axm, color="#6aa9d8", alpha=0.45, lw=0)
            if not mask_r0.is_empty:
                gpd.GeoSeries([mask_r0], crs=holden.crs).plot(
                    ax=axm, color="#1d4e89", alpha=0.9, edgecolor="k", lw=0.4)
            if has_osm:
                holden.osm_polygons.to_crs(holden.crs).boundary.plot(ax=axm, color="#9b1d8a", lw=1.1, ls="--")
            if len(pts):
                pts.plot(ax=axm, color="#1d4e89", markersize=4 if k == 0 else 10, marker="x", lw=0.8)
            dep[~dep["inside_zone"]].plot(ax=axm, color="white", edgecolor="k", markersize=msize, lw=0.7)
            dep[dep["inside_zone"]].plot(ax=axm, color="#d7263d", edgecolor="k", markersize=msize, lw=0.7)
            axm.set_xlim(x0, x1)
            axm.set_ylim(y0, y1)
            axm.set_aspect("equal")
            axm.set_title(title, fontsize=10.5)
            axm.tick_params(labelsize=7)
            axm.ticklabel_format(style="plain", useOffset=False)
            span_km = (x1 - x0) / 1000
            bar = 5000 if span_km > 20 else 1000
            bx, by = x0 + 0.05 * (x1 - x0), y0 + 0.05 * (y1 - y0)
            axm.plot([bx, bx + bar], [by, by], color="k", lw=3)
            axm.text(bx + bar / 2, by + 0.02 * (y1 - y0), f"{bar // 1000} km", ha="center", fontsize=8)
            if k == 0:
                axm.add_patch(plt.Rectangle((windows[1][1][0], windows[1][1][1]), 2 * half, 2 * half,
                                            fill=False, ec="k", lw=1.0, ls="--"))
        handles = [
            Patch(facecolor="#f4a259", edgecolor="#c0561b", label="ASTER strong-anomaly zones"),
            Patch(facecolor="#1d4e89", edgecolor="k", label="USMIN tailings / dump polygons"),
            Patch(facecolor="#6aa9d8", alpha=0.45, label=f"USMIN points + {int(r_hi)} m buffer"),
            Line2D([], [], marker="x", ls="", color="#1d4e89", label="USMIN adit / prospect / pit"),
            Line2D([], [], marker="o", ls="", mfc="#d7263d", mec="k", label="MRDS critical — hit"),
            Line2D([], [], marker="o", ls="", mfc="white", mec="k", label="MRDS critical — miss"),
            Line2D([], [], color="0.3", lw=0.8, label="TIR footprint"),
        ]
        if has_osm:
            handles.insert(3, Line2D([], [], color="#9b1d8a", lw=1.1, ls="--",
                                     label="OSM mine / quarry polygons"))
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8.5,
                   bbox_to_anchor=(0.66, 0.02), frameon=False)
    fig.subplots_adjust(bottom=0.17, top=0.9, left=0.1, right=0.98)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--buffers", type=float, nargs="+", default=[0.0, 100.0],
                    help="Point-feature buffer radii in metres (0 = polygons only)")
    ap.add_argument("--halos", type=float, nargs="*", default=[500.0, 1000.0],
                    help="Dilation radii (m) applied to waste polygons only")
    ap.add_argument("--sites", nargs="+", default=list(FDR_SITES))
    ap.add_argument("--no-osm", action="store_true", help="Skip the Overpass (OSM) layer")
    ap.add_argument("--refresh", action="store_true", help="Re-download cached USMIN / OSM responses")
    ap.add_argument("--skip-figure", action="store_true")
    args = ap.parse_args()
    buffers = sorted(set(args.buffers) | {0.0})
    halos = sorted(set(args.halos))

    rows: list[dict] = []
    holden: SiteInputs | None = None
    for site_id in args.sites:
        print(f"{site_id:22s}", end=" ", flush=True)
        si = load_site_inputs(site_id, REPO_ROOT, CACHE_DIR, fetch_osm=not args.no_osm, refresh=args.refresh)
        if site_id == "holden":
            holden = si
        rows.extend(analyse_site(si, buffers, halos))
        print(f"usmin={si.usmin_status} (pts={rows[-1]['n_usmin_points']}, polys={rows[-1]['n_usmin_polygons']})"
              f"  osm={si.osm_status}  n_crit={len(si.deposits)}")

    df = pd.DataFrame(rows)
    if PHASE4_CSV.is_file():
        ph4 = pd.read_csv(PHASE4_CSV)[["site_id", "n_crit", "hits_crit", "p0_null", "p_binomial"]]
        ph4.columns = ["site_id", "phase4_n_crit", "phase4_hits_crit", "phase4_p0_null", "phase4_p_binomial"]
        df = df.merge(ph4, on="site_id", how="left")
    df["p_binomial"] = df["p_binomial"].astype(float)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {OUT_CSV}  ({len(df)} rows)\n")

    r_hi = int(max(buffers))
    h_cols = [f"halo{int(h)}" for h in halos]
    piv = []
    for s in args.sites:
        d = df[df["site_id"] == s].set_index("variant")
        b = d.loc["baseline"]
        r0 = d.loc["removed_usmin_r0"]
        rh = d.loc[f"removed_usmin_r{r_hi}"]
        x0 = d.loc["excluded_usmin_r0"]
        row = {
            "site": s + (" *" if s in ROBUST_SITES else ""),
            "zone%_r0": r0["pct_zone_on_mask"],
            "fp%_r0": r0["pct_footprint_on_mask"],
            f"zone%_r{r_hi}": rh["pct_zone_on_mask"],
            "E_base": b["enrichment"], "p_base": f"{b['p_binomial']:.2g}",
            "E_rm0": r0["enrichment"], "p_rm0": f"{r0['p_binomial']:.2g}",
        }
        for hc in h_cols:
            hr = d.loc[f"removed_usmin_{hc}"]
            row[f"E_{hc}"] = hr["enrichment"]
            row[f"p_{hc}"] = f"{hr['p_binomial']:.2g}"
        row.update({
            f"E_rm{r_hi}": rh["enrichment"], f"p_rm{r_hi}": f"{rh['p_binomial']:.2g}",
            f"n_rm{r_hi}": int(rh["n_crit"]),
            "p_excl0": f"{x0['p_binomial']:.2g}",
        })
        piv.append(row)
    print(pd.DataFrame(piv).to_string(index=False))
    print("\n(* = fully robust; zone% / fp% = % of zone / footprint area on the USMIN mask;"
          " E = enrichment; rm = mask removed from zones+footprint; haloH = waste polygons dilated H m and"
          " removed; excl = deposits on mask dropped only, p0 unchanged)")

    if "holden" in args.sites and holden is not None:
        d = df[df["site_id"] == "holden"].set_index("variant")
        b, r0, rh = d.loc["baseline"], d.loc["removed_usmin_r0"], d.loc[f"removed_usmin_r{r_hi}"]
        polys, _, _ = mask_layers(holden, "usmin")
        waste = disturbed_mask(polys, None, 0)
        parts = list(getattr(holden.zones_union, "geoms", [holden.zones_union]))
        local = [g for g in parts if g.intersects(waste)]
        in_local = points_in_mask(holden.deposits.geometry, gpd.GeoSeries(local).union_all()) if local else \
            np.zeros(len(holden.deposits), dtype=bool)
        local_hits = int((in_local & holden.deposits["inside_zone"].values).sum())
        largest = max(g.area for g in parts) / holden.zones_union.area * 100
        print("\nHOLDEN (mine-waste hypothesis)")
        print(f"  USMIN polygons in footprint: {int(b['n_usmin_polygons'])} "
              f"({b['usmin_polygon_types'] or 'none'}); points: {int(b['n_usmin_points'])}; OSM: {b['osm_status']}")
        print(f"  Zone area on tailings/dump polygons: {r0['pct_zone_on_mask']:.3f}% "
              f"(footprint on mask {r0['pct_footprint_on_mask']:.3f}%; affinity {r0['zone_mask_affinity']}x)")
        print(f"  Zone patches touching the USMIN waste polygons: {len(local)} "
              f"({sum(g.area for g in local) / holden.zones_union.area * 100:.2f}% of zone area) "
              f"holding {local_hits}/{int(b['hits_crit'])} critical hits; "
              f"largest single zone patch = {largest:.0f}% of zone area")
        print(f"  Critical hits on waste polygons: {int(r0['hits_on_mask'])}/{int(b['hits_crit'])}; "
              f"with {r_hi} m point buffers: {int(rh['hits_on_mask'])}/{int(b['hits_crit'])}")
        worst = [r0["p_binomial"], rh["p_binomial"]]
        line = (f"  Baseline E={b['enrichment']} p={b['p_binomial']:.2g} | polygons removed E={r0['enrichment']} "
                f"p={r0['p_binomial']:.2g}")
        for src in ("usmin", "usmin+osm"):
            for hc in h_cols:
                key = f"removed_{src}_{hc}"
                if key in d.index:
                    hr = d.loc[key]
                    worst.append(hr["p_binomial"])
                    line += f" | {src} {hc} E={hr['enrichment']} p={hr['p_binomial']:.2g}"
        line += f" | +{r_hi} m buffers E={rh['enrichment']} p={rh['p_binomial']:.2g}"
        print(line)
        detected = r0["zone_mask_affinity"] > 1.5
        print("  VERDICT: " + (
            "mine waste plausibly explains the Holden signal." if max(worst) >= 0.05 else
            ("the tailings ARE inside a strong-anomaly patch (ASTER sees the waste), but that local patch holds "
             f"only {local_hits} of {int(b['hits_crit'])} hits; the site-level enrichment survives every mask "
             "(incl. OSM present-day mine outline and 1 km halos), so mine waste does NOT explain Holden's "
             "significance." if detected else
             "mapped mine waste does NOT explain the Holden signal.")))

    lost = []
    for s in sorted(ROBUST_SITES & set(args.sites)):
        d = df[df["site_id"] == s]
        rem = d[d["mode"] == "removed"]
        bad = rem[rem["p_binomial"] >= 0.05]
        if len(bad):
            lost.append(f"{s} ({', '.join(bad['variant'])})")
    print("\nRobust sites losing significance (p>=0.05) after masking: " + ("; ".join(lost) if lost else "none"))

    if not args.skip_figure:
        make_figure(df, holden, buffers, halos, OUT_FIG)
        print(f"Saved {OUT_FIG}")


if __name__ == "__main__":
    main()
