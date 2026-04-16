#!/usr/bin/env python3
"""
04_analysis2_conditional.py — Conditional association testing.

For each of 18 ICVF PGS, fit two models per dataset:

  Unconditional:  logit(IDHmt) ~ ICVF_PGS + age + sex + PC1-PC8
  Conditional:    logit(IDHmt) ~ ICVF_PGS + rs55705857_dosage + age + sex + PC1-PC8

Meta-analyse each model across datasets.
FDR correction across 18 meta-analysis p-values (per model type).
Output both in the same file with a 'model' column.

Usage
-----
python 04_analysis2_conditional.py \
    --data-dir merged_data/ \
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
    RS55705857_COL,
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
        description="Analysis 2: conditional association testing (rs55705857 as covariate)."
    )
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--phenotype", default="idhmt")
    return parser.parse_args()


def main():
    args = parse_args()
    log = setup_logging("analysis2_conditional")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Discover datasets ------------------------------------------------
    datasets = discover_datasets(args.data_dir, args.phenotype)
    if not datasets:
        log.error("No datasets found in %s for phenotype '%s'", args.data_dir, args.phenotype)
        sys.exit(1)
    log.info("Found %d datasets: %s", len(datasets), list(datasets.keys()))

    # --- Define model specifications --------------------------------------
    models = {
        "unconditional": STANDARD_COVARIATES,
        "conditional": [RS55705857_COL] + STANDARD_COVARIATES,
    }

    # --- Per-dataset testing ----------------------------------------------
    rows = []

    for ds_name, ds_path in sorted(datasets.items()):
        df = load_dataset(ds_path)
        log.info("Dataset %s: %d samples", ds_name, len(df))

        for pgs_id in ICVF_PGS_IDS:
            if pgs_id not in df.columns:
                log.warning("Column %s missing in %s — skipping", pgs_id, ds_name)
                continue

            for model_name, covariates in models.items():
                res = run_logistic(df, pgs_id, covariates)
                row = format_result_row(pgs_id, ds_name, model_name, res)
                rows.append(row)

                if not res.get("converged", False):
                    log.warning(
                        "%s / %s / %s: convergence failure (n=%s)",
                        ds_name, pgs_id, model_name, res.get("n"),
                    )

    # --- Meta-analysis per model ------------------------------------------
    for model_name in models:
        for pgs_id in ICVF_PGS_IDS:
            per_ds = [
                r for r in rows
                if r["pgs_id"] == pgs_id
                and r["model"] == model_name
                and r["dataset"] != "meta"
            ]
            meta = meta_analyze_ivw(per_ds)
            mrow = format_meta_row(pgs_id, model_name, meta)
            rows.append(mrow)

    # --- Build output DataFrame -------------------------------------------
    results_df = pd.DataFrame(rows)

    # FDR correction per model type (on meta p-values only)
    for model_name in models:
        mask = (results_df["dataset"] == "meta") & (results_df["model"] == model_name)
        if mask.any():
            results_df.loc[mask, "fdr_q"] = apply_fdr(results_df.loc[mask, "p"].values)

    # Sort
    results_df.sort_values(
        ["model", "pgs_id", "dataset"],
        key=lambda s: s.map({"meta": "zzz_meta"}) if s.name == "dataset" else s,
        inplace=True,
    )
    results_df.reset_index(drop=True, inplace=True)

    # --- Write output -----------------------------------------------------
    out_path = os.path.join(args.outdir, "analysis2_conditional_results.csv")
    results_df.to_csv(out_path, index=False)
    log.info("Results written to %s (%d rows)", out_path, len(results_df))

    # --- Console summary --------------------------------------------------
    print("\n=== Analysis 2: Conditional vs Unconditional (Meta-Analysis) ===\n")
    meta_df = results_df[results_df["dataset"] == "meta"].copy()

    for model_name in ["unconditional", "conditional"]:
        sub = meta_df[meta_df["model"] == model_name]
        if sub.empty:
            continue
        print(f"--- {model_name.upper()} model ---")
        display_cols = ["pgs_id", "beta", "se", "p", "or_", "or_lower", "or_upper",
                        "fdr_q", "I2", "n_datasets"]
        avail_cols = [c for c in display_cols if c in sub.columns]
        print(sub[avail_cols].to_string(index=False))
        print()

    # Side-by-side comparison for the two models
    print("--- COMPARISON: Effect of conditioning on rs55705857 ---")
    uncond = meta_df[meta_df["model"] == "unconditional"].set_index("pgs_id")
    cond = meta_df[meta_df["model"] == "conditional"].set_index("pgs_id")
    common_pgs = uncond.index.intersection(cond.index)
    if len(common_pgs) > 0:
        comparison = pd.DataFrame({
            "pgs_id": common_pgs,
            "beta_uncond": uncond.loc[common_pgs, "beta"].values,
            "p_uncond": uncond.loc[common_pgs, "p"].values,
            "beta_cond": cond.loc[common_pgs, "beta"].values,
            "p_cond": cond.loc[common_pgs, "p"].values,
            "beta_change_pct": (
                (cond.loc[common_pgs, "beta"].values - uncond.loc[common_pgs, "beta"].values)
                / np.where(uncond.loc[common_pgs, "beta"].values != 0,
                           np.abs(uncond.loc[common_pgs, "beta"].values), np.nan)
                * 100
            ),
        })
        print(comparison.to_string(index=False))
        print()

    log.info("Analysis 2 complete.")


if __name__ == "__main__":
    main()
