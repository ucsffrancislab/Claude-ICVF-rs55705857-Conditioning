#!/usr/bin/env python3
"""
05_analysis3_interaction.py — Formal interaction testing.

For each of 18 ICVF PGS, fit per dataset:

  logit(phenotype) ~ ICVF_PGS + rs55705857_dosage
                 + ICVF_PGS × rs55705857_dosage
                 + age + sex + PC1-PC8

Extract:
  - Interaction term coefficient, SE, p-value
  - Main effect of ICVF_PGS from the interaction model

Meta-analyse both the interaction term and the main effect across datasets.
FDR correction across 18 interaction meta-p-values.

Usage
-----
python 05_analysis3_interaction.py \
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
    run_logistic_interaction,
    meta_analyze_ivw,
    apply_fdr,
    format_result_row,
    format_meta_row,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analysis 3: ICVF PGS × rs55705857 interaction testing."
    )
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--phenotype", default="idhmt")
    return parser.parse_args()


def main():
    args = parse_args()
    log = setup_logging("analysis3_interaction")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Discover datasets ------------------------------------------------
    datasets = discover_datasets(args.data_dir, args.phenotype)
    if not datasets:
        log.error("No datasets found in %s for phenotype '%s'", args.data_dir, args.phenotype)
        sys.exit(1)
    log.info("Found %d datasets: %s", len(datasets), list(datasets.keys()))

    # --- Per-dataset interaction testing -----------------------------------
    rows = []

    for ds_name, ds_path in sorted(datasets.items()):
        df = load_dataset(ds_path)
        log.info("Dataset %s: %d samples", ds_name, len(df))

        for pgs_id in ICVF_PGS_IDS:
            if pgs_id not in df.columns:
                log.warning("Column %s missing in %s — skipping", pgs_id, ds_name)
                continue

            res = run_logistic_interaction(
                df, pgs_id, RS55705857_COL, STANDARD_COVARIATES,
            )

            # Main effect of PGS from the interaction model
            main_res = res["main"]
            main_res.update(n=res["n"], n_cases=res["n_cases"],
                            n_controls=res["n_controls"], converged=res["converged"])
            row_main = format_result_row(pgs_id, ds_name, "interaction", main_res, term="pgs_main")
            rows.append(row_main)

            # Interaction term (PGS × rs55705857)
            inter_res = res["interaction_term"]
            inter_res.update(n=res["n"], n_cases=res["n_cases"],
                             n_controls=res["n_controls"], converged=res["converged"])
            row_inter = format_result_row(pgs_id, ds_name, "interaction", inter_res, term="pgs_x_rs55705857")
            rows.append(row_inter)

            # Main effect of rs55705857 (for reference)
            rs_res = res["interaction_var"]
            rs_res.update(n=res["n"], n_cases=res["n_cases"],
                          n_controls=res["n_controls"], converged=res["converged"])
            row_rs = format_result_row(pgs_id, ds_name, "interaction", rs_res, term="rs55705857_main")
            rows.append(row_rs)

            if not res.get("converged", False):
                log.warning(
                    "%s / %s: interaction model convergence failure (n=%s)",
                    ds_name, pgs_id, res.get("n"),
                )

    # --- Meta-analysis per term -------------------------------------------
    for term in ["pgs_main", "pgs_x_rs55705857", "rs55705857_main"]:
        for pgs_id in ICVF_PGS_IDS:
            per_ds = [
                r for r in rows
                if r["pgs_id"] == pgs_id
                and r["term"] == term
                and r["dataset"] != "meta"
            ]
            meta = meta_analyze_ivw(per_ds)
            mrow = format_meta_row(pgs_id, "interaction", meta, term=term)
            rows.append(mrow)

    # --- Build output DataFrame -------------------------------------------
    results_df = pd.DataFrame(rows)

    # FDR correction for interaction term meta-p-values
    mask_inter_meta = (
        (results_df["dataset"] == "meta")
        & (results_df["term"] == "pgs_x_rs55705857")
    )
    if mask_inter_meta.any():
        results_df.loc[mask_inter_meta, "fdr_q"] = apply_fdr(
            results_df.loc[mask_inter_meta, "p"].values
        )

    # Also FDR for main effect
    mask_main_meta = (
        (results_df["dataset"] == "meta")
        & (results_df["term"] == "pgs_main")
    )
    if mask_main_meta.any():
        results_df.loc[mask_main_meta, "fdr_q"] = apply_fdr(
            results_df.loc[mask_main_meta, "p"].values
        )

    # Sort
    results_df.sort_values(
        ["term", "pgs_id", "dataset"],
        key=lambda s: s.map({"meta": "zzz_meta"}) if s.name == "dataset" else s,
        inplace=True,
    )
    results_df.reset_index(drop=True, inplace=True)

    # --- Write output -----------------------------------------------------
    out_path = os.path.join(args.outdir, "analysis3_interaction_results.csv")
    results_df.to_csv(out_path, index=False)
    log.info("Results written to %s (%d rows)", out_path, len(results_df))

    # --- Console summary --------------------------------------------------
    print("\n=== Analysis 3: Interaction Results (Meta-Analysis) ===\n")
    meta_df = results_df[results_df["dataset"] == "meta"].copy()

    print("--- INTERACTION TERM: ICVF_PGS × rs55705857 ---")
    inter_meta = meta_df[meta_df["term"] == "pgs_x_rs55705857"]
    if not inter_meta.empty:
        display_cols = ["pgs_id", "beta", "se", "p", "or_", "or_lower", "or_upper",
                        "fdr_q", "I2", "n_datasets"]
        avail_cols = [c for c in display_cols if c in inter_meta.columns]
        print(inter_meta[avail_cols].to_string(index=False))
    print()

    print("--- MAIN EFFECT: ICVF_PGS (from interaction model) ---")
    main_meta = meta_df[meta_df["term"] == "pgs_main"]
    if not main_meta.empty:
        display_cols = ["pgs_id", "beta", "se", "p", "or_", "or_lower", "or_upper",
                        "fdr_q", "I2", "n_datasets"]
        avail_cols = [c for c in display_cols if c in main_meta.columns]
        print(main_meta[avail_cols].to_string(index=False))
    print()

    log.info("Analysis 3 complete.")


if __name__ == "__main__":
    main()
