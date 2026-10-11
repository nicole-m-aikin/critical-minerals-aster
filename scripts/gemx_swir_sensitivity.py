"""Cutoff sensitivity of the n=2 TIR / SWIR correspondence.

Pre-registered before looking at the grid
-----------------------------------------
Universe: same 90 m overlap as ``gemx_swir_correspondence`` (valid TIR +
SWIR + classification).  No new cubes.  Not Tetracorder.

Primary index: alunite-like band depth.  Kaolinite-like and Al–OH are
reported as secondary.

Cuts: TIR combined score >= {2, 3, 4}; SWIR-high = depth >= overlap
percentile {80, 90, 95}.  The published cell is score >= 3, p90.

Predictions (written before the sweep):

1. Jaccard stays below 0.20 at every nearby cell.
2. No published pin becomes *both* TIR-strong and SWIR-high
   (a coincidence flip).  A pin dropping out of "high" or "strong"
   because the cut got stricter is a threshold dropout, not a fail.
3. At the published cell (3, p90) the already-reported pin story holds:
   Goldfield town = SWIR-high / not TIR-strong;
   Marysvale Alunite Ridge = SWIR-high / not TIR-strong;
   Marysvale town = TIR-strong / not SWIR-high.

If 1–3 fail, the n=2 complementary claim is cutoff-fragile.

Usage:
    conda run -n aster-minerals python scripts/gemx_swir_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from gemx_goldfield_correspondence import (  # noqa: E402
    binary_lift,
    jaccard,
)
from gemx_swir_correspondence import compute  # noqa: E402
from gemx_swir_pilot import SITE_PINS, _site_id_from_argv  # noqa: E402

SCORE_CUTS = (2, 3, 4)
PCT_CUTS = (80, 90, 95)
PUBLISHED = (3, 90)
JACCARD_MAX = 0.20
SITES = ("goldfield", "marysvale")
PRIMARY = "alunite"

PIN_EXPECT_PUBLISHED = {
    "goldfield": {
        "town": {"in_strong_zone": False, "alunite_high": True},
        "target": None,  # district pin is not part of the complementary story
    },
    "marysvale": {
        "town": {"in_strong_zone": True, "alunite_high": False},
        "target": {"in_strong_zone": False, "alunite_high": True},
    },
}


def cutoff_masks(
    overlap: np.ndarray,
    combined: np.ndarray,
    depth: np.ndarray,
    score_min: int,
    pct: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """TIR-strong and SWIR-high masks on the overlap at one cutoff pair."""
    strong = overlap & (combined >= score_min)
    finite = overlap & np.isfinite(depth)
    if not np.any(finite):
        thresh = float("nan")
        high = np.zeros_like(overlap)
    else:
        thresh = float(np.nanpercentile(depth[overlap], pct))
        high = finite & (depth >= thresh)
    return strong, high, thresh


def cutoff_cell(
    overlap: np.ndarray,
    combined: np.ndarray,
    depth: np.ndarray,
    score_min: int,
    pct: float,
) -> dict:
    strong, high, thresh = cutoff_masks(overlap, combined, depth, score_min, pct)
    ev = high[overlap]
    grp = strong[overlap]
    both = int(np.count_nonzero(strong & high))
    return {
        "score_min": int(score_min),
        "swir_pct": int(pct),
        "swir_threshold": thresh,
        "n_overlap": int(np.count_nonzero(overlap)),
        "n_strong": int(np.count_nonzero(strong)),
        "n_high": int(np.count_nonzero(high)),
        "n_both": both,
        "jaccard": float(jaccard(strong, high)),
        "lift": float(binary_lift(ev, grp)),
        "frac_strong": float(np.mean(grp)) if ev.size else None,
    }


def pin_cell(
    overlap: np.ndarray,
    strong: np.ndarray,
    high: np.ndarray,
    depth: np.ndarray,
    row: int,
    col: int,
    label: str,
) -> dict:
    inside = 0 <= row < overlap.shape[0] and 0 <= col < overlap.shape[1]
    depth_v = (
        float(depth[row, col])
        if inside and np.isfinite(depth[row, col])
        else None
    )
    return {
        "label": label,
        "in_overlap": bool(inside and overlap[row, col]),
        "in_strong_zone": bool(inside and strong[row, col]),
        "alunite_high": bool(inside and high[row, col]),
        "alunite": depth_v,
        "both": bool(inside and strong[row, col] and high[row, col]),
    }


def _pin_rc(dst_t, crs: str, lonlat: tuple[float, float], shape: tuple[int, int]):
    from pyproj import Transformer

    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = to_utm.transform(*lonlat)
    c, r = ~dst_t * (x, y)
    r, c = int(r), int(c)
    return r, c


def sweep_site(site) -> dict:
    stats, arrays = compute(site)
    overlap = arrays["overlap"]
    combined = arrays["combined"]
    alun = arrays["alunite"]
    dst_t = arrays["dst_t"]
    crs = arrays["crs"]
    pins_cfg = SITE_PINS[site.id]
    pin_rc = {
        "town": _pin_rc(dst_t, crs, pins_cfg["town"], arrays["shape"]),
        "target": _pin_rc(dst_t, crs, pins_cfg["target"], arrays["shape"]),
    }

    grid = []
    pin_grid = []
    for score_min in SCORE_CUTS:
        for pct in PCT_CUTS:
            cell = cutoff_cell(overlap, combined, alun, score_min, pct)
            strong, high, _ = cutoff_masks(overlap, combined, alun, score_min, pct)
            pins = {}
            for key, (r, c) in pin_rc.items():
                label = pins_cfg["town_label"] if key == "town" else pins_cfg["target_label"]
                pins[key] = pin_cell(overlap, strong, high, alun, r, c, label)
            cell["pins"] = pins
            grid.append(cell)
            pin_grid.append({"score_min": score_min, "swir_pct": pct, "pins": pins})

    return {
        "site_id": site.id,
        "site_name": site.name,
        "n_overlap": stats["n_overlap"],
        "published_score_min": PUBLISHED[0],
        "published_swir_pct": PUBLISHED[1],
        "published_jaccard": stats["indices_in_strong"][PRIMARY]["jaccard_p90"],
        "published_lift": stats["indices_in_strong"][PRIMARY]["lift_p90"],
        "grid": grid,
        "pin_grid": pin_grid,
    }


def evaluate_predictions(site_results: list[dict]) -> dict:
    """Score the pre-registered rules.  Do not move the thresholds after seeing numbers."""
    checks = []
    jaccard_ok = True
    coincidence_ok = True
    published_ok = True

    for site in site_results:
        sid = site["site_id"]
        for cell in site["grid"]:
            jac = cell["jaccard"]
            over = isinstance(jac, float) and jac >= JACCARD_MAX
            if over:
                jaccard_ok = False
            checks.append({
                "rule": "jaccard_below_0.20",
                "site_id": sid,
                "score_min": cell["score_min"],
                "swir_pct": cell["swir_pct"],
                "value": jac,
                "pass": not over,
            })
            for key, pin in cell["pins"].items():
                if pin["both"]:
                    coincidence_ok = False
                    checks.append({
                        "rule": "no_coincidence_flip",
                        "site_id": sid,
                        "pin": key,
                        "score_min": cell["score_min"],
                        "swir_pct": cell["swir_pct"],
                        "pass": False,
                    })

        expect = PIN_EXPECT_PUBLISHED[sid]
        published_cell = next(
            c for c in site["grid"]
            if c["score_min"] == PUBLISHED[0] and c["swir_pct"] == PUBLISHED[1]
        )
        for key, want in expect.items():
            if want is None:
                continue
            got = published_cell["pins"][key]
            ok = (
                got["in_strong_zone"] == want["in_strong_zone"]
                and got["alunite_high"] == want["alunite_high"]
            )
            if not ok:
                published_ok = False
            checks.append({
                "rule": "published_pin_story",
                "site_id": sid,
                "pin": key,
                "expected": want,
                "got": {
                    "in_strong_zone": got["in_strong_zone"],
                    "alunite_high": got["alunite_high"],
                },
                "pass": ok,
            })

    return {
        "jaccard_stays_below_0.20": jaccard_ok,
        "no_coincidence_flip": coincidence_ok,
        "published_pin_story_holds": published_ok,
        "cutoff_fragile": not (jaccard_ok and coincidence_ok and published_ok),
        "checks": checks,
        "note": (
            "Pre-registered.  Complementary means the maps do not occupy the "
            "same pixels, not that SWIR explains the epithermal TIR nulls."
        ),
    }


def _grid_matrix(site: dict, key: str) -> np.ndarray:
    mat = np.full((len(SCORE_CUTS), len(PCT_CUTS)), np.nan)
    for cell in site["grid"]:
        i = SCORE_CUTS.index(cell["score_min"])
        j = PCT_CUTS.index(cell["swir_pct"])
        mat[i, j] = cell[key]
    return mat


def save_figure(site_results: list[dict], verdict: dict) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.2), constrained_layout=True)
    if not verdict["jaccard_stays_below_0.20"] or not verdict["no_coincidence_flip"]:
        tag = "cutoff-fragile (Jaccard or coincidence)"
    elif not verdict["published_pin_story_holds"]:
        tag = "Jaccard holds nearby · Goldfield town pin is not SWIR-high"
    else:
        tag = "complementary holds nearby"
    title = f"TIR / SWIR cutoff sensitivity  ·  n=2  ·  not mineral ID  ·  {tag}"
    fig.suptitle(title, fontsize=12)

    for row, site in enumerate(site_results):
        jac = _grid_matrix(site, "jaccard")
        lift = _grid_matrix(site, "lift")
        for col, (mat, label, vmax, cmap) in enumerate((
            (jac, "Jaccard", 0.25, "YlOrRd"),
            (lift, "Lift", 2.5, "RdYlBu_r"),
        )):
            ax = axes[row, col]
            im = ax.imshow(mat, vmin=0.0, vmax=vmax, cmap=cmap, origin="upper")
            ax.set_xticks(range(len(PCT_CUTS)), [f"p{p}" for p in PCT_CUTS])
            ax.set_yticks(range(len(SCORE_CUTS)), [f"score ≥ {s}" for s in SCORE_CUTS])
            ax.set_title(f"{site['site_id']}  {label}", fontsize=10)
            for i in range(len(SCORE_CUTS)):
                for j in range(len(PCT_CUTS)):
                    val = mat[i, j]
                    if not np.isfinite(val):
                        txt = "—"
                    elif label == "Jaccard":
                        txt = f"{val:.3f}"
                    else:
                        txt = f"{val:.2f}×"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                            color="black")
            pi = SCORE_CUTS.index(PUBLISHED[0])
            pj = PCT_CUTS.index(PUBLISHED[1])
            ax.add_patch(Rectangle(
                (pj - 0.5, pi - 0.5), 1, 1, fill=False,
                edgecolor="black", linewidth=1.8,
            ))
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    axes[0, 0].set_ylabel("TIR cutoff")
    axes[1, 0].set_ylabel("TIR cutoff")
    axes[1, 0].set_xlabel("SWIR percentile on the overlap")
    axes[1, 1].set_xlabel("SWIR percentile on the overlap")
    fig.text(
        0.5, -0.02,
        "Black box = published cell (score ≥ 3, p90).  "
        f"Fail if any Jaccard ≥ {JACCARD_MAX:.2f} or a pin is both TIR-strong and SWIR-high.",
        ha="center", fontsize=8, color="#333333",
    )
    out = ROOT / "figures" / "gemx_swir_sensitivity.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main(site_ids: tuple[str, ...] = SITES) -> dict:
    from critical_minerals_aster.config import load_site_by_id

    site_results = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, ROOT / "sites")
        print(f"sweep {site_id} …", flush=True)
        site_results.append(sweep_site(site))
    verdict = evaluate_predictions(site_results)
    payload = {
        "hypothesis": (
            "complementary (not coincident) is not an artifact of "
            "score ≥ 3 and overlap p90"
        ),
        "score_cuts": list(SCORE_CUTS),
        "swir_percentiles": list(PCT_CUTS),
        "jaccard_max": JACCARD_MAX,
        "sites": site_results,
        "verdict": {
            k: v for k, v in verdict.items() if k != "checks"
        },
        "checks": verdict["checks"],
        "note": (
            "Reuse of existing Goldfield and Marysvale GeoTIFFs.  "
            "Not Tetracorder.  Not a 24-site test."
        ),
    }
    stats_path = ROOT / "results" / "gemx_swir_sensitivity.json"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(payload, indent=2))
    fig_path = save_figure(site_results, verdict)
    print(json.dumps(payload["verdict"], indent=2))
    print(f"Wrote {stats_path.relative_to(ROOT)}")
    print(f"Wrote {fig_path.relative_to(ROOT)}")
    return payload


if __name__ == "__main__":
    ids = SITES
    if "--site" in sys.argv:
        ids = (_site_id_from_argv(),)
    main(ids)
