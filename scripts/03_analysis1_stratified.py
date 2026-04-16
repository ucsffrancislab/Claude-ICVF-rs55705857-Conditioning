#!/usr/bin/env python3
"""
03_analysis1_stratified.py — Stratified association testing.

Split each dataset by rs55705857 carrier status:
  - Non-carriers (AA genotype): dosage < 0.5
  - Carriers: dosage >= 0.5

Within each stratum, test each ICVF PGS:
  logit(IDHmt) ~ ICVF_PGS + age + sex + PC1-PC8

Meta-analyse across datasets within each stratum.
FDR correction across the 18 meta-analysis p-values (per stratum).

Usage
-----
python 03_analysis1_stratified.py \
    --data-dir merged_data/ \
    --outdir results/ \
    --phenotype idhmt
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

# Allow running from the project root
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
        description="Analysis 1: stratified association testing by rs55705857 carrier status."
    )
    parser.add_argument(
        "--data-dir", required=True,
        help="Directory containing per-dataset analysis-ready TSV files.",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory for output files.",
    )
    parser.add_argument(
        "--phenotype", default="idhmt",
        help="Phenotype label used in file naming (default: idhmt).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    log = setup_logging("analysis1_stratified")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Discover datasets ------------------------------------------------
    datasets = discover_datasets(args.data_dir, args.phenotype)
    if not datasets:
        log.error("No datasets found in %s for phenotype '%s'", args.data_dir, args.phenotype)
        sys.exit(1)
    log.info("Found %d datasets: %s", len(datasets), list(datasets.keys()))

    # --- Stratification & testing -----------------------------------------
    rows = []
    strata = {"noncarrier": "dosage < 0.5", "carrier": "dosage >= 0.5"}

    # Collect sample-size summary
    sample_summary = []

    for ds_name, ds_path in sorted(datasets.items()):
        df = load_dataset(ds_path)
        log.info("Dataset %s: %d samples", ds_name, len(df))

        # Stratify
        mask_noncarrier = df[RS55705857_COL] < 0.5
        mask_carrier = df[RS55705857_COL] >= 0.5
        strata_data = {
            "noncarrier": df.loc[mask_noncarrier].copy(),
            "carrier": df.loc[mask_carrier].copy(),
        }

        for stratum_name, sub_df in strata_data.items():
            n_total = len(sub_df)
            # Handle both PLINK (2/1) and standard (1/0) phenotype coding
            if sub_df["phenotype"].max() > 1:
                n_cases = int((sub_df["phenotype"] == 2).sum())
                n_controls = int((sub_df["phenotype"] == 1).sum())
            else:
                n_cases = int((sub_df["phenotype"] == 1).sum())
                n_controls = int((sub_df["phenotype"] == 0).sum())
            sample_summary.append(dict(
                dataset=ds_name, stratum=stratum_name,
                n_total=n_total, n_cases=n_cases, n_controls=n_controls,
            ))
            log.info(
                "  %s / %s: n=%d (cases=%d, controls=%d)",
                ds_name, stratum_name, n_total, n_cases, n_controls,
            )

            for pgs_id in ICVF_PGS_IDS:
                if pgs_id not in sub_df.columns:
                    log.warning("Column %s missing in %s — skipping", pgs_id, ds_name)
                    continue

                res = run_logistic(sub_df, pgs_id, STANDARD_COVARIATES)
                row = format_result_row(pgs_id, ds_name, f"stratified_{stratum_name}", res)
                row["stratum"] = stratum_name
                rows.append(row)

    # --- Meta-analysis per stratum ----------------------------------------
    for stratum_name in ["noncarrier", "carrier"]:
        model_label = f"stratified_{stratum_name}"
        for pgs_id in ICVF_PGS_IDS:
            per_ds = [
                r for r in rows
                if r["pgs_id"] == pgs_id
                and r["model"] == model_label
                and r["dataset"] != "meta"
            ]
            meta = meta_analyze_ivw(per_ds)
            mrow = format_meta_row(pgs_id, model_label, meta)
            mrow["stratum"] = stratum_name
            rows.append(mrow)

    # --- Build output DataFrame -------------------------------------------
    results_df = pd.DataFrame(rows)

    # FDR correction per stratum (on meta-analysis p-values only)
    for stratum_name in ["noncarrier", "carrier"]:
        mask = (results_df["dataset"] == "meta") & (results_df["stratum"] == stratum_name)
        if mask.any():
            results_df.loc[mask, "fdr_q"] = apply_fdr(results_df.loc[mask, "p"].values)

    # Sort for readability
    results_df.sort_values(
        ["stratum", "pgs_id", "dataset"],
        key=lambda s: s.map({"meta": "zzz_meta"}) if s.name == "dataset" else s,
        inplace=True,
    )
    results_df.reset_index(drop=True, inplace=True)

    # --- Write outputs ----------------------------------------------------
    out_results = os.path.join(args.outdir, "analysis1_stratified_results.csv")
    results_df.to_csv(out_results, index=False)
    log.info("Results written to %s (%d rows)", out_results, len(results_df))

    # Sample size summary
    out_samples = os.path.join(args.outdir, "analysis1_stratified_sample_sizes.csv")
    pd.DataFrame(sample_summary).to_csv(out_samples, index=False)
    log.info("Sample size summary written to %s", out_samples)

    # --- Console summary --------------------------------------------------
    print("\n=== Analysis 1: Stratified Association Results (Meta-Analysis) ===\n")
    meta_df = results_df[results_df["dataset"] == "meta"].copy()
    for stratum_name in ["noncarrier", "carrier"]:
        sub = meta_df[meta_df["stratum"] == stratum_name]
        if sub.empty:
            continue
        print(f"--- {stratum_name.upper()} stratum ---")
        display_cols = ["pgs_id", "beta", "se", "p", "or_", "or_lower", "or_upper",
                        "fdr_q", "I2", "n_datasets"]
        avail_cols = [c for c in display_cols if c in sub.columns]
        print(sub[avail_cols].to_string(index=False))
        print()

    log.info("Analysis 1 complete.")


if __name__ == "__main__":
    main()
