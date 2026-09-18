"""Aggregate per-site summaries into a national comparison table."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_site_summaries(
    results_dir: Path,
    row_types: list[str] | None = None,
) -> pd.DataFrame:
    """Load all *_summary.csv files, optionally filtered by row_type.

    Parameters
    ----------
    row_types:
        If given, keep only rows whose ``row_type`` is in this list.
        Defaults to ``["site"]`` to preserve the original behaviour for
        callers that expect one row per site.  Pass ``None`` to return all
        rows (site + commodity + earth_mri).
    """
    _ALL = object()  # sentinel: include every row_type
    _filter = _ALL if row_types is None else row_types
    results_dir = Path(results_dir)
    frames: list[pd.DataFrame] = []
    for path in sorted(
        p for p in results_dir.glob("*_summary.csv") if "national" not in p.name
    ):
        df = pd.read_csv(path)
        if "row_type" in df.columns and _filter is not _ALL:
            df = df[df["row_type"].isin(_filter)]  # type: ignore[arg-type]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


#: Colours for the two FDR-corrected significance tiers plus the null tier.
SIG_PRIMARY_COLOR = "#c0392b"    # red   — BH-FDR significant, exact binomial test
SIG_SECONDARY_COLOR = "#e67e22"  # orange — BH-FDR significant on Mann-Whitney U only
SIG_NONE_COLOR = "#5d8aa8"       # blue-grey — not significant under FDR


def load_fdr_significance(results_dir: Path) -> pd.DataFrame | None:
    """Load FDR-corrected per-site significance, indexed by ``site_id``.

    Reads ``results/phase4_fdr_corrected_significance.csv`` — the output of the
    corrected framework described in docs/results.md — and exposes the two
    significance tiers the synthesis figures colour by:

    ``sig_primary``
        Benjamini-Hochberg-corrected one-sided exact binomial test against the
        site's own geometric null.  This is the headline claim.
    ``sig_secondary_only``
        Clears BH-FDR on the threshold-free Mann-Whitney U test on the
        continuous 0-6 score, but *not* on the binomial test.

    Returns ``None`` when the file is absent, so callers can fall back to
    all-deposit hit rates without significance colouring.
    """
    path = Path(results_dir) / "phase4_fdr_corrected_significance.csv"
    if not path.exists():
        return None
    fdr = pd.read_csv(path)
    binomial = fdr["sig_fdr_binomial"].fillna(False).astype(bool)
    mannwhitney = fdr["sig_fdr_mannwhitney"].fillna(False).astype(bool)
    out = pd.DataFrame(
        {
            "site_id": fdr["site_id"],
            "hr_crit_pct": fdr["hit_rate_crit_pct"],
            "n_crit": fdr["n_crit"],
            "sig_primary": binomial,
            "sig_secondary_only": mannwhitney & ~binomial,
        }
    )
    return out.set_index("site_id")


def _tier_color(row) -> str:
    if row["sig_primary"]:
        return SIG_PRIMARY_COLOR
    if row["sig_secondary_only"]:
        return SIG_SECONDARY_COLOR
    return SIG_NONE_COLOR


def _per_site_barh_geometry(n: int) -> tuple[tuple[float, float], float]:
    """Figure size and y-tick font size for a one-row-per-site horizontal chart.

    Per-row height shrinks once the survey passes ~50 sites so the figure keeps
    a roughly 2:1 aspect ratio.  Without this the 109-site run rendered a
    10 x 65 inch strip (6.6:1) whose title and legend were unreadable.
    """
    row_h = 0.6 if n <= 50 else max(0.2, 30.0 / n)
    height = max(4.0, n * row_h)
    width = 10.0 if n <= 50 else 13.0
    fontsize = 8.0 if n <= 50 else 6.0
    return (width, height), fontsize


def _fdr_legend_handles(n_primary: int, n_secondary: int):
    from matplotlib.patches import Patch

    return [
        Patch(
            facecolor=SIG_PRIMARY_COLOR, edgecolor="black", linewidth=0.5,
            label=f"FDR-significant, binomial (n={n_primary})",
        ),
        Patch(
            facecolor=SIG_SECONDARY_COLOR, edgecolor="black", linewidth=0.5,
            label=f"FDR-significant, continuous score only (n={n_secondary})",
        ),
        Patch(
            facecolor=SIG_NONE_COLOR, edgecolor="black", linewidth=0.5,
            label="Not significant under FDR",
        ),
    ]


def save_paired_hitrate_figure(results_dir: Path, figures_dir: Path) -> Path | None:
    """Generate figures/07_hitrate_comparison.png — all-deposit vs critical-only hit rate per site.

    Each site gets two dots connected by a line: grey for the all-deposit rate
    (potentially inflated by aggregate/stone) and coloured for the critical-mineral
    rate, tiered by FDR-corrected significance.  Sites are sorted by all-deposit
    rate so the inflation story is visible as leftward shifts.
    """
    results_dir = Path(results_dir)
    fdr = load_fdr_significance(results_dir)
    if fdr is None:
        return None

    national = load_site_summaries(results_dir, row_types=["site"])
    if national.empty:
        return None

    sig = (
        national.set_index("site_id")[["site_name", "hit_rate_pct"]]
        .rename(columns={"hit_rate_pct": "hr_all_pct"})
        .join(fdr, how="inner")
        .dropna(subset=["hr_crit_pct"])
        .sort_values("hr_all_pct", ascending=True)
        .reset_index()
    )
    if sig.empty:
        return None

    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    n = len(sig)
    y = np.arange(n)
    figsize, tick_fontsize = _per_site_barh_geometry(n)

    fig, ax = plt.subplots(figsize=figsize)

    # Connecting lines
    for i, row in sig.iterrows():
        ax.plot(
            [row["hr_all_pct"], row["hr_crit_pct"]],
            [i, i],
            color="#dddddd",
            linewidth=1.2,
            zorder=1,
        )

    # All-deposits dots (grey, open)
    ax.scatter(
        sig["hr_all_pct"], y,
        color="#999999", s=55, zorder=2, label="All deposits",
        edgecolors="#555555", linewidths=0.6,
    )

    crit_colors = [_tier_color(row) for _, row in sig.iterrows()]
    ax.scatter(
        sig["hr_crit_pct"], y,
        c=crit_colors, s=65, zorder=3, label="Critical minerals only",
        edgecolors="black", linewidths=0.5,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(sig["site_id"], fontsize=tick_fontsize)
    ax.set_ylim(-1, n)
    ax.set_xlabel("Hit rate (% of deposits in strong TIR zones)")
    n_primary = int(sig["sig_primary"].sum())
    n_secondary = int(sig["sig_secondary_only"].sum())
    ax.set_title(
        "All-deposit vs critical-mineral hit rates by site\n"
        f"{n} sites  ·  grey = all deposits, coloured = critical minerals only  ·  "
        f"significance is BH-FDR corrected (α = 0.05)"
    )
    ax.axvline(0, color="black", linewidth=0.5, alpha=0.4)
    ax.grid(axis="x", alpha=0.25, linewidth=0.5)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#999999", edgecolor="#555555", linewidth=0.6, label="All deposits"),
        *_fdr_legend_handles(n_primary, n_secondary),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="lower right")

    plt.tight_layout()
    out = figures_dir / "07_hitrate_comparison.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")
    return out


def write_national_summary(results_dir: Path) -> Path:
    """Write national_summary.csv, national_summary.parquet, and results.duckdb under results/."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    national = load_site_summaries(results_dir, row_types=None)
    csv_path = results_dir / "national_summary.csv"
    national.to_csv(csv_path, index=False)
    parquet_path = results_dir / "national_summary.parquet"
    try:
        national.to_parquet(parquet_path, index=False)
    except ImportError:
        parquet_path = None

    # --- Task 3: DuckDB output ---
    try:
        import duckdb
        db_path = results_dir / "results.duckdb"
        con = duckdb.connect(str(db_path))
        con.execute("DROP TABLE IF EXISTS site_summaries")
        con.execute("CREATE TABLE site_summaries AS SELECT * FROM national")
        con.close()
        print(f"  DuckDB: {db_path}")
    except ImportError:
        pass
    except Exception as exc:
        print(f"  DuckDB write failed: {exc}", file=sys.stderr)

    figures_dir = results_dir.parent / "figures"

    # --- Task 1: structure-hit rate scatter ---
    save_structure_hitrate_scatter(national, results_dir, figures_dir)

    # --- Task 2: national stacked bar (Earth MRI categories, critical only) ---
    save_national_figure(results_dir, figures_dir)

    # --- Task 3: HTML figure index ---
    write_figure_index(results_dir, figures_dir)

    # --- Task 4: paired all-deposits vs critical-only hit rate comparison ---
    save_paired_hitrate_figure(results_dir, figures_dir)

    return csv_path


def save_structure_hitrate_scatter(
    national_df: pd.DataFrame, results_dir: Path, figures_dir: Path
) -> None:
    """Generate figures/06_structure_hit_rate.png — structural proximity vs critical-mineral hit rate."""
    figures_dir = Path(figures_dir)
    results_dir = Path(results_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    df = national_df[national_df["row_type"] == "site"].copy()
    df = df.dropna(subset=["mean_nearest_structure_m"])

    if len(df) < 2:
        print(
            "Warning: fewer than 2 sites have structure data; "
            "skipping figure 06_structure_hit_rate.png"
        )
        return

    # Merge critical-mineral hit rates and FDR-corrected significance tiers.
    # Sites absent from the FDR table are dropped rather than defaulted to
    # "not significant" — silently recolouring untested sites grey previously
    # mislabelled 64 of 109 sites when the significance table lagged the survey.
    fdr = load_fdr_significance(results_dir)
    use_crit = False
    if fdr is not None:
        df = df.set_index("site_id").join(fdr, how="inner").reset_index()
        df = df.dropna(subset=["hr_crit_pct"])
        if len(df) < 2:
            print(
                "Warning: fewer than 2 sites have both structure and FDR data; "
                "skipping figure 06_structure_hit_rate.png"
            )
            return
        df["hr_plot"] = df["hr_crit_pct"]
        use_crit = True
    else:
        df["hr_plot"] = df["hit_rate_pct"]
        df["sig_primary"] = False
        df["sig_secondary_only"] = False

    out = figures_dir / "06_structure_hit_rate.png"

    sizes = np.clip(
        df["n_deposits_bbox"] / df["n_deposits_bbox"].max() * 400, 30, 400
    )

    point_colors = [_tier_color(row) for _, row in df.iterrows()]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(
        df["mean_nearest_structure_m"] / 1000,
        df["hr_plot"],
        s=sizes,
        c=point_colors,
        alpha=0.85,
        edgecolors="black",
        linewidths=0.5,
        zorder=3,
    )

    # Label only the FDR-significant sites once the survey is large enough that
    # annotating every point turns the dense mid-field into unreadable overlap.
    if use_crit and len(df) > 60:
        labelled = df[df["sig_primary"] | df["sig_secondary_only"]]
    else:
        labelled = df
    for _, row in labelled.iterrows():
        ax.annotate(
            row["site_id"],
            (row["mean_nearest_structure_m"] / 1000, row["hr_plot"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    mean_hr = df["hr_plot"].mean()
    ax.axhline(
        mean_hr,
        color="#2c3e50",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label=f"Mean hit rate ({mean_hr:.1f}%)",
    )

    n_primary = int(df["sig_primary"].sum())
    n_secondary = int(df["sig_secondary_only"].sum())
    legend_elements = [ax.lines[0], *_fdr_legend_handles(n_primary, n_secondary)]
    ax.legend(handles=legend_elements, fontsize=8)

    ax.set_xscale("log")
    y_label = (
        "Critical-mineral hit rate (% of critical deposits in strong TIR zones)"
        if use_crit else "MRDS deposit hit rate (%)"
    )
    ax.set_xlabel("Mean distance to nearest fault (km, log scale)")
    ax.set_ylabel(y_label)
    label_note = (
        "  ·  FDR-significant sites labelled" if len(labelled) < len(df) else ""
    )
    ax.set_title(
        f"Structural proximity vs spectral detectability\n"
        f"{len(df)} ASTER study sites with fault data  ·  "
        f"{n_primary} FDR-significant (binomial)  ·  {n_secondary} continuous-score only"
        f"{label_note}"
    )
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


def write_figure_index(results_dir: Path, figures_dir: Path) -> Path:
    """Generate figures/index.html — a sortable grid of site cards with key metrics."""
    results_dir = Path(results_dir)
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    csv_path = results_dir / "national_summary.csv"
    df = pd.read_csv(csv_path)
    df = df[df["row_type"] == "site"].sort_values("hit_rate_pct", ascending=False)
    n_sites = len(df)

    cards_html = []
    for _, row in df.iterrows():
        site_id = row["site_id"]
        site_name = row.get("site_name", site_id)
        hit_rate = row.get("hit_rate_pct", 0.0)
        n_in = row.get("n_deposits_in_zones", 0)
        n_bbox = row.get("n_deposits_bbox", 0)
        n_zones = row.get("n_zones", 0)

        struct_html = ""
        mean_struct = row.get("mean_nearest_structure_m")
        n_on_struct = row.get("n_deposits_on_structure")
        try:
            if not pd.isna(mean_struct) and not pd.isna(n_on_struct):
                struct_html = (
                    f'<p class="struct">{int(n_on_struct)} on structure '
                    f"({mean_struct / 1000:.1f} km mean)</p>"
                )
        except (TypeError, ValueError):
            pass

        card = f"""    <div class="card" data-hit-rate="{hit_rate}" data-deposits="{n_bbox}" data-name="{site_id}">
      <img src="sites/{site_id}/03_deposit_overlay.png" alt="{site_name}" loading="lazy">
      <div class="card-body">
        <h3>{site_name}</h3>
        <p class="stats">Hit rate: {hit_rate:.1f}% &middot; {int(n_in)}/{int(n_bbox)} deposits &middot; {int(n_zones)} zones</p>
        {struct_html}
      </div>
    </div>"""
        cards_html.append(card)

    generated_ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    cards_joined = "\n".join(cards_html)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Critical Minerals ASTER &mdash; Site Overview</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: #f4f6f9;
      color: #333;
      padding: 1.5rem;
    }}
    header {{
      max-width: 1200px;
      margin: 0 auto 1.5rem;
    }}
    h1 {{ font-size: 1.6rem; margin-bottom: 0.5rem; color: #1a1a2e; }}
    .sort-controls {{
      display: flex;
      gap: 0.5rem;
      flex-wrap: wrap;
      margin-bottom: 1.25rem;
    }}
    .sort-controls button {{
      padding: 0.4rem 0.9rem;
      border: 1px solid #ccc;
      border-radius: 4px;
      background: #fff;
      cursor: pointer;
      font-size: 0.85rem;
      transition: background 0.15s, border-color 0.15s;
    }}
    .sort-controls button:hover, .sort-controls button.active {{
      background: #2c3e50;
      color: #fff;
      border-color: #2c3e50;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
      gap: 1rem;
      max-width: 1200px;
      margin: 0 auto;
    }}
    .card {{
      background: #fff;
      border-radius: 8px;
      overflow: hidden;
      box-shadow: 0 1px 4px rgba(0,0,0,0.1);
      transition: transform 0.15s, box-shadow 0.15s;
    }}
    .card:hover {{
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }}
    .card img {{
      width: 100%;
      display: block;
      aspect-ratio: 4/3;
      object-fit: cover;
      background: #e8ecf0;
    }}
    .card-body {{ padding: 0.75rem 1rem 1rem; }}
    .card-body h3 {{ font-size: 0.95rem; margin-bottom: 0.35rem; color: #1a1a2e; }}
    .card-body .stats {{ font-size: 0.78rem; color: #555; margin-bottom: 0.25rem; }}
    .card-body .struct {{ font-size: 0.75rem; color: #777; }}
    footer {{
      text-align: center;
      margin-top: 2rem;
      font-size: 0.78rem;
      color: #888;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Critical Minerals ASTER &mdash; Site Overview</h1>
    <div class="sort-controls">
      <button class="active" onclick="sortCards('hit-rate')">Sort by hit rate</button>
      <button onclick="sortCards('deposits')">Sort by deposits</button>
      <button onclick="sortCards('name')">Sort by site name</button>
    </div>
  </header>
  <div class="grid" id="grid">
{cards_joined}
  </div>
  <footer>Generated {generated_ts} &middot; {n_sites} sites</footer>
  <script>
    function sortCards(key) {{
      const grid = document.getElementById('grid');
      const cards = Array.from(grid.querySelectorAll('.card'));
      cards.sort((a, b) => {{
        if (key === 'name') return a.dataset.name.localeCompare(b.dataset.name);
        if (key === 'deposits') return Number(b.dataset.deposits) - Number(a.dataset.deposits);
        return Number(b.dataset.hitRate) - Number(a.dataset.hitRate);
      }});
      cards.forEach(c => grid.appendChild(c));
      document.querySelectorAll('.sort-controls button').forEach(btn => btn.classList.remove('active'));
      event.target.classList.add('active');
    }}
  </script>
</body>
</html>
"""

    out = figures_dir / "index.html"
    out.write_text(html, encoding="utf-8")
    print(f"  Figure index: {out}")
    return out


def save_national_figure(results_dir: Path, figures_dir: Path) -> Path:
    """Generate figures/05_national_hit_rates.png — stacked bar chart by Earth MRI category.

    Bars are broken down by Earth MRI commodity category so the figure shows
    both hit rate and deposit composition in one view.  Uses Paul Tol's Muted
    palette (colorblind-safe for up to 10 categories).
    """
    results_dir = Path(results_dir)
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    national = load_site_summaries(results_dir, row_types=["site"])
    out = figures_dir / "05_national_hit_rates.png"

    if national.empty:
        return out

    # Load earth_mri rows, excluding Non-Critical so bars reflect critical-mineral signal only
    earth_mri = load_site_summaries(results_dir, row_types=["earth_mri"])
    earth_mri = earth_mri[earth_mri["earth_mri_category"] != "Non-Critical"]

    # Use critical-mineral hit rates for bar length and ordering when available
    fdr = load_fdr_significance(results_dir)
    if fdr is not None:
        site_id_to_name = national.set_index("site_id")["site_name"]
        crit_by_name = (
            fdr.assign(site_name=site_id_to_name.reindex(fdr.index))
            .dropna(subset=["site_name"])
            .set_index("site_name")["hr_crit_pct"]
        )
        name_order = national.set_index("site_name").index
        crit_hr = crit_by_name.reindex(name_order).fillna(
            national.set_index("site_name")["hit_rate_pct"]
        )
        site_order = crit_hr.sort_values(ascending=True).index.tolist()
        site_hit_rate = crit_hr
    else:
        site_order = national.sort_values("hit_rate_pct", ascending=True)["site_name"].tolist()
        site_hit_rate = national.set_index("site_name")["hit_rate_pct"]

    pivot = earth_mri.pivot_table(
        index="site_name", columns="earth_mri_category",
        values="n_deposits_in_zones", aggfunc="sum", fill_value=0,
    )
    pivot = pivot.reindex(site_order).fillna(0)

    row_totals = pivot.sum(axis=1).replace(0, 1)
    cat_share = pivot.div(row_totals, axis=0)
    hr_aligned = site_hit_rate.reindex(pivot.index)
    pivot_abs = cat_share.mul(hr_aligned, axis=0)

    # Paul Tol "Muted" palette — colorblind-safe for up to 10 categories.
    EARTH_MRI_COLORS = {
        "Gold/Silver":                  "#DDCC77",
        "Base Metals":                  "#88CCEE",
        "Non-Critical":                 "#DDDDDD",
        "Energy":                       "#CC6677",
        "Battery Metals \u2013 Co/Ni":  "#44AA99",
        "Battery Metals \u2013 Li/Brine": "#117733",
        "Specialty/High-Tech":          "#332288",
        "Industrial":                   "#999933",
        "REE":                          "#AA4499",
        "PGM":                          "#882255",
    }
    CAT_ORDER = [
        "Gold/Silver", "Base Metals", "Non-Critical", "Energy",
        "Battery Metals \u2013 Co/Ni", "Battery Metals \u2013 Li/Brine",
        "Specialty/High-Tech", "Industrial", "REE", "PGM",
    ]
    cat_cols = [c for c in CAT_ORDER if c in pivot_abs.columns]
    cat_cols += [c for c in pivot_abs.columns if c not in cat_cols]

    figsize, tick_fontsize = _per_site_barh_geometry(len(site_order))
    fig, ax = plt.subplots(figsize=figsize)
    left = pd.Series(0.0, index=pivot_abs.index)
    for cat in cat_cols:
        vals = pivot_abs[cat]
        ax.barh(pivot_abs.index, vals, left=left,
                color=EARTH_MRI_COLORS.get(cat, "#AAAAAA"),
                label=cat, edgecolor="white", linewidth=0.4)
        left = left + vals

    x_max = max(hr_aligned.max() * 1.08, 22)
    ax.set_xlim(0, x_max)
    ax.tick_params(axis="y", labelsize=tick_fontsize)
    ax.set_xlabel("Critical-mineral hit rate (% of critical deposits in strong TIR zones)")
    ax.set_title(
        f"Alteration\u2013deposit correlation by site, {len(site_order)} sites "
        "(critical minerals only)\n"
        "(Non-Critical excluded  \u00b7  colour = Earth MRI category of in-zone critical deposits)"
    )
    ax.legend(
        loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0,
        fontsize=7.5, title="Earth MRI category", title_fontsize=8,
        framealpha=0.9,
    )
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out
