#!/usr/bin/env python3
"""
07_figures.py — Publication-quality figures for ICVF conditioning analysis.

Generates 4 figures (PNG + PDF, 300 DPI):
  1. Forest plot: conditional vs unconditional ORs (Analysis 2)
  2. Effect comparison scatter: IDHwt OR vs IDHmt conditional OR
  3. Stratified forest plot: AA non-carrier stratum (Analysis 1)
  4. LD pruning comparison: original vs pruned effect sizes

Usage
-----
python 07_figures.py \
    --results-dir results/ \
    --outdir figures/ \
    --idh-wt-results results/idhwt_results.csv   # optional
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PGS_SHORT_NAMES = {
    "PGS001454": "Body CC",
    "PGS001456": "Cerebral peduncle R",
    "PGS001457": "Cing. cing. gyrus L",
    "PGS001458": "Cing. cing. gyrus R",
    "PGS001459": "Cing. hippocampus L",
    "PGS001460": "Cing. hippocampus R",
    "PGS001466": "Genu CC",
    "PGS001471": "Mid. cerebellar ped.",
    "PGS001474": "Post. int. capsule L",
    "PGS001478": "Retrolent. int. cap. L",
    "PGS001479": "Retrolent. int. cap. R",
    "PGS001480": "Sagittal stratum L",
    "PGS001481": "Sagittal stratum R",
    "PGS001484": "Sup. corona rad. L",
    "PGS001485": "Sup. corona rad. R",
    "PGS001662": "Acoustic rad. R\u2020",
    "PGS001669": "Forceps major\u2020",
    "PGS001679": "Parahipp. cing. R\u2020",
}

NEGATIVE_CONTROLS = {"PGS001456", "PGS001662", "PGS001679"}

# Consistent color scheme
COLOR_UNCOND = "#2166ac"      # blue
COLOR_COND = "#b2182b"        # red
COLOR_PRUNED = "#4daf4a"      # green
COLOR_STRAT = "#984ea3"       # purple
COLOR_NEG_CTRL = "#ff7f00"    # orange for negative controls
COLOR_FDR_SIG = "#d62728"     # red filled
COLOR_NS = "#7f7f7f"          # grey for non-significant
COLOR_IDHWT = "#1b9e77"       # teal

DPI = 300
SINGLE_COL = 3.5   # inches
DOUBLE_COL = 7.0   # inches

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate publication figures for ICVF conditioning analysis."
    )
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--idh-wt-results", default=None,
                        help="CSV with IDHwt unconditional results for scatter plot.")
    return parser.parse_args()


def short_name(pgs_id):
    return PGS_SHORT_NAMES.get(pgs_id, pgs_id)


def load_meta_results(filepath, model_filter=None, term_filter=None):
    """Load a results CSV and return only meta-analysis rows, optionally filtered."""
    df = pd.read_csv(filepath)
    mask = df["dataset"] == "meta"
    if model_filter is not None:
        mask &= df["model"] == model_filter
    if term_filter is not None:
        mask &= df["term"] == term_filter
    return df.loc[mask].copy().reset_index(drop=True)


def save_fig(fig, outdir, name):
    """Save figure as both PNG and PDF."""
    fig.savefig(os.path.join(outdir, f"{name}.png"))
    fig.savefig(os.path.join(outdir, f"{name}.pdf"))
    plt.close(fig)
    print(f"  Saved {name}.png and {name}.pdf")



# =========================================================================
# Figure 1: Forest plot — conditional vs unconditional
# =========================================================================

def figure1_forest_conditional(results_dir, outdir):
    """Two-panel forest plot: unconditional (left) vs conditional (right)."""
    filepath = os.path.join(results_dir, "analysis2_conditional_results.csv")
    if not os.path.isfile(filepath):
        print("  WARNING: analysis2_conditional_results.csv not found — skipping Fig 1")
        return

    uncond = load_meta_results(filepath, model_filter="unconditional")
    cond = load_meta_results(filepath, model_filter="conditional")

    if uncond.empty or cond.empty:
        print("  WARNING: No meta-analysis rows found — skipping Fig 1")
        return

    # Sort by conditional p-value
    cond = cond.sort_values("p").reset_index(drop=True)
    sort_order = cond["pgs_id"].tolist()
    uncond = uncond.set_index("pgs_id").loc[sort_order].reset_index()

    n_pgs = len(sort_order)
    y_pos = np.arange(n_pgs)

    fig, (ax_u, ax_c) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 0.32 * n_pgs + 1.2),
        sharey=True, gridspec_kw={"wspace": 0.08},
    )

    for ax, df, color, title in [
        (ax_u, uncond, COLOR_UNCOND, "Unconditional"),
        (ax_c, cond, COLOR_COND, "Conditional on rs55705857"),
    ]:
        for i, (_, row) in enumerate(df.iterrows()):
            pgs_id = row["pgs_id"]
            or_val = row.get("or_", np.nan)
            or_lo = row.get("or_lower", np.nan)
            or_hi = row.get("or_upper", np.nan)
            fdr_q = row.get("fdr_q", 1.0)
            is_neg_ctrl = pgs_id in NEGATIVE_CONTROLS

            if not np.isfinite(or_val):
                continue

            marker_color = COLOR_NEG_CTRL if is_neg_ctrl else color
            marker = "D" if is_neg_ctrl else "o"
            fill = marker_color if (fdr_q < 0.05 and np.isfinite(fdr_q)) else "white"
            edge = marker_color

            ax.plot(or_val, i, marker=marker, color=fill,
                    markeredgecolor=edge, markeredgewidth=1.0, markersize=6, zorder=3)
            ax.plot([or_lo, or_hi], [i, i], color=edge, linewidth=1.0, zorder=2)

        ax.axvline(x=1, color="grey", linestyle="--", linewidth=0.7, zorder=1)
        ax.set_xlabel("Odds Ratio (95% CI)")
        ax.set_title(title, fontweight="bold")

        # Set symmetric x-limits
        all_lo = df["or_lower"].dropna()
        all_hi = df["or_upper"].dropna()
        if not all_lo.empty and not all_hi.empty:
            x_margin = 0.05
            x_lo = max(0.3, all_lo.min() - x_margin)
            x_hi = min(3.0, all_hi.max() + x_margin)
            ax.set_xlim(x_lo, x_hi)

    # Y-axis labels on left panel only
    labels = [short_name(pid) for pid in sort_order]
    ax_u.set_yticks(y_pos)
    ax_u.set_yticklabels(labels)
    ax_u.invert_yaxis()

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_UNCOND,
               markeredgecolor=COLOR_UNCOND, markersize=6, label="FDR < 0.05 (filled)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor="grey", markersize=6, label="Not significant"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_NEG_CTRL,
               markeredgecolor=COLOR_NEG_CTRL, markersize=6,
               label="Negative control (\u2020)"),
    ]
    ax_c.legend(handles=legend_elements, loc="lower right", framealpha=0.9, fontsize=6)

    save_fig(fig, outdir, "fig1_forest_conditional_vs_unconditional")



# =========================================================================
# Figure 2: Effect comparison scatter (IDHwt OR vs IDHmt conditional OR)
# =========================================================================

def figure2_scatter_comparison(results_dir, outdir, idh_wt_path):
    """Scatter: IDHwt OR (x) vs IDHmt conditional OR (y)."""
    if idh_wt_path is None or not os.path.isfile(idh_wt_path):
        print("  INFO: IDHwt results not provided — skipping Fig 2")
        return

    cond_path = os.path.join(results_dir, "analysis2_conditional_results.csv")
    if not os.path.isfile(cond_path):
        print("  WARNING: analysis2_conditional_results.csv not found — skipping Fig 2")
        return

    cond = load_meta_results(cond_path, model_filter="conditional")

    # Load IDHwt results — expect meta-analysis rows with pgs_id, or_, fdr_q
    wt_df = pd.read_csv(idh_wt_path)
    wt_meta = wt_df[wt_df["dataset"] == "meta"].copy()
    if "model" in wt_meta.columns:
        wt_meta = wt_meta[wt_meta["model"] == "unconditional"]

    if cond.empty or wt_meta.empty:
        print("  WARNING: Insufficient data for scatter — skipping Fig 2")
        return

    # Merge on pgs_id
    merged = cond.set_index("pgs_id")[["or_", "or_lower", "or_upper", "fdr_q"]].rename(
        columns={"or_": "or_mt", "or_lower": "or_mt_lo", "or_upper": "or_mt_hi", "fdr_q": "fdr_mt"}
    ).join(
        wt_meta.set_index("pgs_id")[["or_", "or_lower", "or_upper", "fdr_q"]].rename(
            columns={"or_": "or_wt", "or_lower": "or_wt_lo", "or_upper": "or_wt_hi", "fdr_q": "fdr_wt"}
        ),
        how="inner",
    )

    if merged.empty:
        print("  WARNING: No overlapping PGS IDs — skipping Fig 2")
        return

    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))

    for pgs_id, row in merged.iterrows():
        is_neg = pgs_id in NEGATIVE_CONTROLS
        sig_mt = row.get("fdr_mt", 1) < 0.05
        sig_wt = row.get("fdr_wt", 1) < 0.05

        if is_neg:
            c, marker = COLOR_NEG_CTRL, "D"
        elif sig_mt or sig_wt:
            c, marker = COLOR_FDR_SIG, "o"
        else:
            c, marker = COLOR_NS, "o"

        ax.plot(row["or_wt"], row["or_mt"], marker=marker, color=c,
                markeredgecolor=c, markersize=6, zorder=3)

        # Error bars
        ax.plot([row["or_wt_lo"], row["or_wt_hi"]], [row["or_mt"], row["or_mt"]],
                color=c, linewidth=0.5, alpha=0.4, zorder=2)
        ax.plot([row["or_wt"], row["or_wt"]], [row["or_mt_lo"], row["or_mt_hi"]],
                color=c, linewidth=0.5, alpha=0.4, zorder=2)

        # Label notable points (most extreme)
        if abs(row["or_mt"] - 1) > 0.15 or abs(row["or_wt"] - 1) > 0.15 or is_neg:
            ax.annotate(short_name(pgs_id), (row["or_wt"], row["or_mt"]),
                        fontsize=5, ha="left", va="bottom",
                        xytext=(3, 3), textcoords="offset points")

    # Identity line
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "k--", linewidth=0.7, alpha=0.5, zorder=1)

    # Reference lines at OR=1
    ax.axhline(y=1, color="grey", linewidth=0.5, alpha=0.3)
    ax.axvline(x=1, color="grey", linewidth=0.5, alpha=0.3)

    # Pearson r
    valid = merged.dropna(subset=["or_wt", "or_mt"])
    if len(valid) > 2:
        from scipy.stats import pearsonr
        r, r_p = pearsonr(valid["or_wt"], valid["or_mt"])
        ax.text(0.05, 0.95, f"r = {r:.2f} (p = {r_p:.3f})",
                transform=ax.transAxes, fontsize=7, va="top")

    ax.set_xlabel("IDHwt OR (unconditional)")
    ax.set_ylabel("IDHmt OR (conditional on rs55705857)")
    ax.set_title("Effect Comparison: IDHwt vs IDHmt", fontweight="bold")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_FDR_SIG,
               markeredgecolor=COLOR_FDR_SIG, markersize=5, label="FDR < 0.05"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_NS,
               markeredgecolor=COLOR_NS, markersize=5, label="Not significant"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_NEG_CTRL,
               markeredgecolor=COLOR_NEG_CTRL, markersize=5, label="Negative control"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=6)

    save_fig(fig, outdir, "fig2_scatter_idhwt_vs_idhmt")



# =========================================================================
# Figure 3: Stratified forest plot (AA non-carriers)
# =========================================================================

def figure3_forest_stratified(results_dir, outdir):
    """Single-panel forest: ORs in the AA non-carrier stratum."""
    filepath = os.path.join(results_dir, "analysis1_stratified_results.csv")
    if not os.path.isfile(filepath):
        print("  WARNING: analysis1_stratified_results.csv not found — skipping Fig 3")
        return

    df_all = pd.read_csv(filepath)
    aa = df_all[(df_all["dataset"] == "meta") & (df_all["stratum"] == "noncarrier")].copy()

    if aa.empty:
        print("  WARNING: No AA non-carrier meta rows — skipping Fig 3")
        return

    aa = aa.sort_values("p").reset_index(drop=True)
    n_pgs = len(aa)
    y_pos = np.arange(n_pgs)

    fig, ax = plt.subplots(figsize=(SINGLE_COL, 0.32 * n_pgs + 1.0))

    for i, (_, row) in enumerate(aa.iterrows()):
        pgs_id = row["pgs_id"]
        or_val = row.get("or_", np.nan)
        or_lo = row.get("or_lower", np.nan)
        or_hi = row.get("or_upper", np.nan)
        fdr_q = row.get("fdr_q", 1.0)
        is_neg = pgs_id in NEGATIVE_CONTROLS

        if not np.isfinite(or_val):
            continue

        mc = COLOR_NEG_CTRL if is_neg else COLOR_STRAT
        marker = "D" if is_neg else "o"
        fill = mc if (fdr_q < 0.05 and np.isfinite(fdr_q)) else "white"

        ax.plot(or_val, i, marker=marker, color=fill,
                markeredgecolor=mc, markeredgewidth=1.0, markersize=6, zorder=3)
        ax.plot([or_lo, or_hi], [i, i], color=mc, linewidth=1.0, zorder=2)

    ax.axvline(x=1, color="grey", linestyle="--", linewidth=0.7, zorder=1)
    labels = [short_name(pid) for pid in aa["pgs_id"]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Odds Ratio (95% CI)")
    ax.set_title("AA Non-Carrier Stratum", fontweight="bold")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_STRAT,
               markeredgecolor=COLOR_STRAT, markersize=6, label="FDR < 0.05"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor=COLOR_STRAT, markersize=6, label="Not significant"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_NEG_CTRL,
               markeredgecolor=COLOR_NEG_CTRL, markersize=6, label="Negative control"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", framealpha=0.9, fontsize=6)

    save_fig(fig, outdir, "fig3_forest_stratified_noncarrier")


# =========================================================================
# Figure 4: LD pruning comparison
# =========================================================================

def figure4_ld_pruning(results_dir, outdir):
    """Paired dot plot: unconditional vs LD-pruned effect sizes."""
    uncond_path = os.path.join(results_dir, "analysis2_conditional_results.csv")
    pruned_path = os.path.join(results_dir, "analysis4_ld_pruned_results.csv")

    if not os.path.isfile(uncond_path) or not os.path.isfile(pruned_path):
        print("  WARNING: Missing analysis2 or analysis4 results — skipping Fig 4")
        return

    uncond = load_meta_results(uncond_path, model_filter="unconditional")
    pruned = load_meta_results(pruned_path, model_filter="ld_pruned")

    if uncond.empty or pruned.empty:
        print("  WARNING: Empty meta rows — skipping Fig 4")
        return

    # Merge
    merged = uncond.set_index("pgs_id")[["or_", "or_lower", "or_upper"]].rename(
        columns={"or_": "or_orig", "or_lower": "lo_orig", "or_upper": "hi_orig"}
    ).join(
        pruned.set_index("pgs_id")[["or_", "or_lower", "or_upper"]].rename(
            columns={"or_": "or_pruned", "or_lower": "lo_pruned", "or_upper": "hi_pruned"}
        ),
        how="inner",
    )

    if merged.empty:
        print("  WARNING: No overlapping PGS — skipping Fig 4")
        return

    # Sort by original OR
    merged = merged.sort_values("or_orig")
    n_pgs = len(merged)
    y_pos = np.arange(n_pgs)
    offset = 0.15

    fig, ax = plt.subplots(figsize=(SINGLE_COL, 0.36 * n_pgs + 1.2))

    for i, (pgs_id, row) in enumerate(merged.iterrows()):
        is_neg = pgs_id in NEGATIVE_CONTROLS

        # Original (top of pair)
        mc_o = COLOR_NEG_CTRL if is_neg else COLOR_UNCOND
        ax.plot(row["or_orig"], i - offset, marker="o", color=mc_o,
                markeredgecolor=mc_o, markersize=5, zorder=3)
        ax.plot([row["lo_orig"], row["hi_orig"]], [i - offset, i - offset],
                color=mc_o, linewidth=0.8, zorder=2)

        # Pruned (bottom of pair)
        mc_p = COLOR_NEG_CTRL if is_neg else COLOR_PRUNED
        marker_p = "D" if is_neg else "s"
        ax.plot(row["or_pruned"], i + offset, marker=marker_p, color=mc_p,
                markeredgecolor=mc_p, markersize=5, zorder=3)
        ax.plot([row["lo_pruned"], row["hi_pruned"]], [i + offset, i + offset],
                color=mc_p, linewidth=0.8, zorder=2)

        # Connecting line
        ax.plot([row["or_orig"], row["or_pruned"]], [i - offset, i + offset],
                color="grey", linewidth=0.4, alpha=0.6, zorder=1)

    ax.axvline(x=1, color="grey", linestyle="--", linewidth=0.7, zorder=0)

    labels = [short_name(pid) for pid in merged.index]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Odds Ratio (95% CI)")
    ax.set_title("Original vs LD-Pruned PGS", fontweight="bold")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_UNCOND,
               markeredgecolor=COLOR_UNCOND, markersize=5, label="Original"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=COLOR_PRUNED,
               markeredgecolor=COLOR_PRUNED, markersize=5, label="LD-pruned"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_NEG_CTRL,
               markeredgecolor=COLOR_NEG_CTRL, markersize=5, label="Negative control"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", framealpha=0.9, fontsize=6)

    save_fig(fig, outdir, "fig4_ld_pruning_comparison")



# =========================================================================
# Main
# =========================================================================

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("=== Generating publication figures ===")
    print(f"  Results dir: {args.results_dir}")
    print(f"  Output dir:  {args.outdir}")
    print()

    print("[Figure 1] Forest plot: conditional vs unconditional")
    figure1_forest_conditional(args.results_dir, args.outdir)

    print("[Figure 2] Scatter: IDHwt vs IDHmt conditional")
    figure2_scatter_comparison(args.results_dir, args.outdir, args.idh_wt_results)

    print("[Figure 3] Forest plot: AA non-carrier stratum")
    figure3_forest_stratified(args.results_dir, args.outdir)

    print("[Figure 4] LD pruning comparison")
    figure4_ld_pruning(args.results_dir, args.outdir)

    print()
    print("=== All figures complete ===")


if __name__ == "__main__":
    main()
