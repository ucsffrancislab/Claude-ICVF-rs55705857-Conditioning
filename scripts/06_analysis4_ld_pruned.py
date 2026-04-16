#!/usr/bin/env python3
"""
06_analysis4_ld_pruned.py — LD-pruned PGS association testing.

Read pre-computed LD-pruned ICVF PGS scores (from 06a_prep_ld_pruning)
and run unconditional association tests:

  logit(IDHmt) ~ ICVF_PGS_pruned + age + sex + PC1-PC8

Meta-analyse across datasets.  FDR correction across 18 tests.
Compare to Analysis 2 unconditional results if available.

The pruned-scores directory should contain per-dataset files with
columns: IID plus one column per pruned PGS (named {PGS_ID}_pruned).

Usage
-----
python 06_analysis4_ld_pruned.py \
    --data-dir merged_data/ \
    --pruned-scores-dir ld_pruned_scores/ \
    --outdir results/ \
    --phenotype idhmt
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.analysis_utils import (
    STANDARD_COVARIATES,
    ICVF_PGS_IDS,
    PGS_WITHOUT_RS55705857,
    setup_logging,
    load_dataset,
    discover_datasets,
    run_logistic,
    meta_analyze_ivw,
    apply_fdr,
    format_result_row,
    format_meta_row,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analysis 4: LD-pruned PGS association testing."
    )
    parser.add_argument("--data-dir", required=True,
                        help="Directory with per-dataset analysis-ready TSVs.")
    parser.add_argument("--pruned-scores-dir", required=True,
                        help="Directory with LD-pruned PGS score files from 06a.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--phenotype", default="idhmt")
    parser.add_argument("--analysis2-results", default=None,
                        help="Path to analysis2_conditional_results.csv for comparison.")
    return parser.parse_args()


def load_pruned_scores(pruned_dir: str, dataset_name: str) -> pd.DataFrame | None:
    """Load pruned PGS scores for a dataset.

    Expected file: {dataset_name}_pruned_scores.tsv
    Columns: IID, {PGS_ID}_pruned (z-scored pruned PGS values)
    """
    # Try several naming conventions
    candidates = [
        f"{dataset_name}_pruned_scores.tsv",
        f"{dataset_name}_pruned_pgs.tsv",
        f"{dataset_name}_ld_pruned_scores.tsv",
    ]
    for fname in candidates:
        fpath = os.path.join(pruned_dir, fname)
        if os.path.isfile(fpath):
            return pd.read_csv(fpath, sep="\t")
    return None


def main():
    args = parse_args()
    log = setup_logging("analysis4_ld_pruned")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Discover datasets ------------------------------------------------
    datasets = discover_datasets(args.data_dir, args.phenotype)
    if not datasets:
        log.error("No datasets found in %s for phenotype '%s'", args.data_dir, args.phenotype)
        sys.exit(1)
    log.info("Found %d datasets: %s", len(datasets), list(datasets.keys()))

    # --- Load pruning summary (for reporting) -----------------------------
    pruning_summary_path = os.path.join(args.pruned_scores_dir, "pruning_summary.tsv")
    pruning_summary = None
    if os.path.isfile(pruning_summary_path):
        pruning_summary = pd.read_csv(pruning_summary_path, sep="\t")
        log.info("Loaded pruning summary: %d PGS", len(pruning_summary))

    # --- Per-dataset testing with pruned scores ---------------------------
    rows = []

    for ds_name, ds_path in sorted(datasets.items()):
        df_base = load_dataset(ds_path)
        log.info("Dataset %s: %d samples", ds_name, len(df_base))

        # Load pruned scores
        pruned_df = load_pruned_scores(args.pruned_scores_dir, ds_name)
        if pruned_df is None:
            log.warning("No pruned scores found for %s — skipping", ds_name)
            continue

        # Merge on IID
        df = df_base.merge(pruned_df, on="IID", how="inner")
        log.info("  After merging pruned scores: %d samples", len(df))

        for pgs_id in ICVF_PGS_IDS:
            pruned_col = f"{pgs_id}_pruned"

            # For PGS without rs55705857, the pruned score equals the original
            # Still test for consistency, but note this in the output
            if pruned_col not in df.columns:
                # Fall back to original if no pruning was needed
                if pgs_id in PGS_WITHOUT_RS55705857 and pgs_id in df.columns:
                    log.info(
                        "  %s: rs55705857 not in scoring file; using original PGS",
                        pgs_id,
                    )
                    df[pruned_col] = df[pgs_id]
                else:
                    log.warning("Column %s not found for %s — skipping", pruned_col, ds_name)
                    continue

            res = run_logistic(df, pruned_col, STANDARD_COVARIATES)
            row = format_result_row(pgs_id, ds_name, "ld_pruned", res)
            row["pgs_col_used"] = pruned_col
            row["rs55705857_in_scoring"] = pgs_id not in PGS_WITHOUT_RS55705857
            rows.append(row)

    # --- Meta-analysis ----------------------------------------------------
    for pgs_id in ICVF_PGS_IDS:
        per_ds = [
            r for r in rows
            if r["pgs_id"] == pgs_id and r["dataset"] != "meta"
        ]
        meta = meta_analyze_ivw(per_ds)
        mrow = format_meta_row(pgs_id, "ld_pruned", meta)
        mrow["rs55705857_in_scoring"] = pgs_id not in PGS_WITHOUT_RS55705857
        rows.append(mrow)

    # --- Build output DataFrame -------------------------------------------
    results_df = pd.DataFrame(rows)

    # FDR correction on meta p-values
    mask_meta = results_df["dataset"] == "meta"
    if mask_meta.any():
        results_df.loc[mask_meta, "fdr_q"] = apply_fdr(
            results_df.loc[mask_meta, "p"].values
        )

    # Sort
    results_df.sort_values(
        ["pgs_id", "dataset"],
        key=lambda s: s.map({"meta": "zzz_meta"}) if s.name == "dataset" else s,
        inplace=True,
    )
    results_df.reset_index(drop=True, inplace=True)

    # --- Write output -----------------------------------------------------
    out_path = os.path.join(args.outdir, "analysis4_ld_pruned_results.csv")
    results_df.to_csv(out_path, index=False)
    log.info("Results written to %s (%d rows)", out_path, len(results_df))

    # --- Optional comparison to Analysis 2 --------------------------------
    if args.analysis2_results and os.path.isfile(args.analysis2_results):
        log.info("Loading Analysis 2 results for comparison: %s", args.analysis2_results)
        a2 = pd.read_csv(args.analysis2_results)
        a2_uncond_meta = a2[
            (a2["dataset"] == "meta") & (a2["model"] == "unconditional")
        ].set_index("pgs_id")

        pruned_meta = results_df[results_df["dataset"] == "meta"].set_index("pgs_id")
        common = pruned_meta.index.intersection(a2_uncond_meta.index)

        if len(common) > 0:
            comparison = pd.DataFrame({
                "pgs_id": common,
                "beta_original": a2_uncond_meta.loc[common, "beta"].values,
                "p_original": a2_uncond_meta.loc[common, "p"].values,
                "beta_pruned": pruned_meta.loc[common, "beta"].values,
                "p_pruned": pruned_meta.loc[common, "p"].values,
                "rs55705857_in_scoring": [
                    pgs not in PGS_WITHOUT_RS55705857 for pgs in common
                ],
            })
            comp_path = os.path.join(args.outdir, "analysis4_vs_analysis2_comparison.csv")
            comparison.to_csv(comp_path, index=False)
            log.info("Comparison written to %s", comp_path)

            print("\n--- Comparison: LD-Pruned vs Original (unconditional) ---")
            print(comparison.to_string(index=False))
            print()

    # --- Console summary --------------------------------------------------
    print("\n=== Analysis 4: LD-Pruned PGS Results (Meta-Analysis) ===\n")
    meta_df = results_df[results_df["dataset"] == "meta"].copy()
    display_cols = ["pgs_id", "beta", "se", "p", "or_", "or_lower", "or_upper",
                    "fdr_q", "I2", "n_datasets", "rs55705857_in_scoring"]
    avail_cols = [c for c in display_cols if c in meta_df.columns]
    print(meta_df[avail_cols].to_string(index=False))
    print()

    # Highlight negative controls
    neg_controls = meta_df[meta_df.get("rs55705857_in_scoring", pd.Series(dtype=bool)) == False]
    if not neg_controls.empty:
        print("--- NEGATIVE CONTROLS (rs55705857 not in scoring file) ---")
        print("These PGS should show minimal change vs original scores:")
        print(neg_controls[["pgs_id", "beta", "p", "fdr_q"]].to_string(index=False))
        print()

    # Attach pruning summary if available
    if pruning_summary is not None:
        print("--- PRUNING SUMMARY ---")
        print(pruning_summary.to_string(index=False))
        print()

    log.info("Analysis 4 complete.")


if __name__ == "__main__":
    main()
