#!/usr/bin/env python3
"""Aggregate per-site results into results/national_summary.csv and regenerate figure 05."""

from pathlib import Path

from critical_minerals_aster.synthesis import save_national_figure, write_national_summary

if __name__ == "__main__":
    repo = Path(__file__).resolve().parents[1]
    out = write_national_summary(repo / "results")
    print(f"Wrote {out}")
    fig_out = save_national_figure(repo / "results", repo / "figures")
    print(f"Saved {fig_out}")
