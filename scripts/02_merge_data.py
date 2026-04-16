#!/usr/bin/env python3
"""
02_merge_data.py

Merge covariates, PGS z-scores, and rs55705857 dosage into analysis-ready
files for the ICVF conditioning analysis.

For each dataset (cidr, i370, onco, tcga):
  1. Load covariates CSV  → drop exclude==1
  2. Load PGS z-scores    → keep only the 18 ICVF-relevant columns
  3. Load rs55705857 dosage TSV
  4. Inner-merge on IID
  5. Define phenotype (idhmt / idhmt_intact / idhmt_codel)
  6. Encode sex as numeric (F=0, M=1)
  7. Drop samples with any missing covariates, PGS, or dosage
  8. Standardize each PGS to mean=0, sd=1 within the dataset
  9. Write {dataset}_analysis_ready.tsv

Also writes a summary TSV with per-dataset sample counts.

Usage:
    python 02_merge_data.py \
        --scores-dir /path/to/scores \
        --dosage-dir /path/to/dosage \
        --outdir      /path/to/out \
        --phenotype   idhmt
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Constants ────────────────────────────────────────────────────────────────

DATASETS = ["cidr", "i370", "onco", "tcga"]

# The 18 ICVF-relevant PGS model IDs
ICVF_PGS_COLUMNS = [
    "PGS001454", "PGS001456", "PGS001457", "PGS001458", "PGS001459",
    "PGS001460", "PGS001466", "PGS001471", "PGS001474", "PGS001478",
    "PGS001479", "PGS001480", "PGS001481", "PGS001484", "PGS001485",
    "PGS001662", "PGS001669", "PGS001679",
]

# Covariates required to be non-missing for inclusion
REQUIRED_COVARIATES = ["age", "sex", "PC1", "PC2", "PC3", "PC4",
                       "PC5", "PC6", "PC7", "PC8"]

DOSAGE_COLUMN = "rs55705857_dosage"

# Valid phenotype definitions
VALID_PHENOTYPES = {"idhmt", "idhmt_intact", "idhmt_codel"}

# ── Logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Phenotype definition ────────────────────────────────────────────────────


def define_phenotype(df: pd.DataFrame, phenotype: str) -> pd.DataFrame:
    """
    Create a binary 'phenotype' column and drop ambiguous / irrelevant samples.

    Parameters
    ----------
    df : DataFrame with at least columns 'case', 'idh', and optionally 'pq'.
    phenotype : one of 'idhmt', 'idhmt_intact', 'idhmt_codel'.

    Returns
    -------
    DataFrame with a new 'phenotype' column (1 = case, 0 = control) and only
    the samples that are unambiguously in the case or control group.
    """
    out = df.copy()

    # Controls are always case==0 (population controls)
    is_control = out["case"] == 0

    if phenotype == "idhmt":
        # Cases: IDH-mutant glioma  (idh == 1)
        is_case = out["idh"] == 1
    elif phenotype == "idhmt_intact":
        # Cases: IDH-mutant, 1p/19q intact  (idh == 1 & pq == 0)
        is_case = (out["idh"] == 1) & (out["pq"] == 0)
    elif phenotype == "idhmt_codel":
        # Cases: IDH-mutant, 1p/19q co-deleted  (idh == 1 & pq == 1)
        is_case = (out["idh"] == 1) & (out["pq"] == 1)
    else:
        raise ValueError(f"Unknown phenotype: {phenotype!r}. "
                         f"Choose from {VALID_PHENOTYPES}")

    # Keep only unambiguous cases and controls
    keep = is_case | is_control
    out = out.loc[keep].copy()
    out["phenotype"] = np.where(is_case.loc[keep], 1, 0)

    return out


# ── Per-dataset processing ──────────────────────────────────────────────────


def process_dataset(
    dataset: str,
    scores_dir: Path,
    dosage_dir: Path,
    outdir: Path,
    phenotype: str,
) -> dict:
    """
    Process one dataset: load, merge, filter, standardize, save.

    Returns a summary dict with sample counts.
    """
    summary = {
        "dataset": dataset,
        "n_covariates_loaded": 0,
        "n_after_exclude": 0,
        "n_pgs_loaded": 0,
        "n_dosage_loaded": 0,
        "n_merged": 0,
        "n_after_phenotype": 0,
        "n_after_missing_drop": 0,
        "n_cases": 0,
        "n_controls": 0,
        "status": "OK",
    }

    log.info("=" * 60)
    log.info(f"Processing dataset: {dataset}")
    log.info("=" * 60)

    # ── 1. Load covariates ───────────────────────────────────────────────
    cov_path = scores_dir / f"{dataset}-covariates.csv"
    if not cov_path.exists():
        log.error(f"Covariate file not found: {cov_path}")
        summary["status"] = "MISSING_COVARIATES"
        return summary

    cov = pd.read_csv(cov_path)
    # Normalise the sample-ID column name to IID
    if "IID" not in cov.columns and "iid" in [c.lower() for c in cov.columns]:
        id_col = [c for c in cov.columns if c.lower() == "iid"][0]
        cov = cov.rename(columns={id_col: "IID"})
    cov["IID"] = cov["IID"].astype(str)
    summary["n_covariates_loaded"] = len(cov)
    log.info(f"  Covariates loaded: {len(cov)} samples")

    # Drop excluded samples
    if "exclude" in cov.columns:
        n_before = len(cov)
        cov = cov[cov["exclude"] != 1].copy()
        n_excluded = n_before - len(cov)
        log.info(f"  Dropped {n_excluded} excluded samples (exclude==1)")
    summary["n_after_exclude"] = len(cov)

    # ── 2. Load PGS z-scores (only the 18 ICVF columns) ─────────────────
    pgs_path = scores_dir / f"{dataset}.scores.z-scores.txt.gz"
    if not pgs_path.exists():
        log.error(f"PGS scores not found: {pgs_path}")
        summary["status"] = "MISSING_PGS"
        return summary

    # Read just the header to find the ID column and target PGS columns
    pgs_header = pd.read_csv(pgs_path, sep=",", nrows=0)
    pgs_cols = list(pgs_header.columns)

    # Identify the sample-ID column (first column; may be "sample" or "IID")
    id_col_pgs = pgs_cols[0]

    # Determine which of the 18 ICVF columns are present
    available_pgs = [c for c in ICVF_PGS_COLUMNS if c in pgs_cols]
    missing_pgs = [c for c in ICVF_PGS_COLUMNS if c not in pgs_cols]
    if missing_pgs:
        log.warning(f"  Missing PGS columns in {dataset}: {missing_pgs}")
    if not available_pgs:
        log.error(f"  No ICVF PGS columns found in {dataset}")
        summary["status"] = "NO_PGS_COLUMNS"
        return summary

    # Read only the columns we need
    usecols = [id_col_pgs] + available_pgs
    pgs = pd.read_csv(pgs_path, sep=",", usecols=usecols)
    pgs = pgs.rename(columns={id_col_pgs: "IID"})
    pgs["IID"] = pgs["IID"].astype(str)
    summary["n_pgs_loaded"] = len(pgs)
    log.info(f"  PGS z-scores loaded: {len(pgs)} samples, "
             f"{len(available_pgs)}/{len(ICVF_PGS_COLUMNS)} PGS columns")

    # ── 3. Load rs55705857 dosage ────────────────────────────────────────
    dos_path = dosage_dir / f"{dataset}_rs55705857_dosage.tsv"
    if not dos_path.exists():
        log.error(f"Dosage file not found: {dos_path}")
        summary["status"] = "MISSING_DOSAGE"
        return summary

    dos = pd.read_csv(dos_path, sep="\t")
    if "IID" not in dos.columns:
        # Try to find an ID column
        dos = dos.rename(columns={dos.columns[0]: "IID"})
    dos["IID"] = dos["IID"].astype(str)
    summary["n_dosage_loaded"] = len(dos)
    log.info(f"  Dosage loaded: {len(dos)} samples")

    # ── 4. Merge on IID ──────────────────────────────────────────────────
    merged = cov.merge(pgs, on="IID", how="inner")
    n_cov_pgs = len(merged)
    merged = merged.merge(dos, on="IID", how="inner")
    summary["n_merged"] = len(merged)
    log.info(f"  After merge: {len(merged)} samples "
             f"(cov∩pgs={n_cov_pgs}, cov∩pgs∩dos={len(merged)})")

    if len(merged) == 0:
        log.error(f"  No samples remain after merge for {dataset}")
        summary["status"] = "EMPTY_AFTER_MERGE"
        return summary

    # ── 5. Define phenotype ──────────────────────────────────────────────
    merged = define_phenotype(merged, phenotype)
    summary["n_after_phenotype"] = len(merged)
    log.info(f"  After phenotype filter ({phenotype}): {len(merged)} samples")

    if len(merged) == 0:
        log.error(f"  No samples remain after phenotype definition for {dataset}")
        summary["status"] = "EMPTY_AFTER_PHENOTYPE"
        return summary

    # ── 6. Encode sex as numeric (F=0, M=1) ─────────────────────────────
    if "sex" in merged.columns:
        sex_map = {"F": 0, "M": 1, "f": 0, "m": 1}
        merged["sex"] = merged["sex"].map(sex_map)
        n_sex_missing = merged["sex"].isna().sum()
        if n_sex_missing > 0:
            log.warning(f"  {n_sex_missing} samples have unmapped sex values")

    # ── 7. Drop samples with missing covariates, PGS, or dosage ─────────
    check_cols = (
        REQUIRED_COVARIATES
        + available_pgs
        + [DOSAGE_COLUMN]
    )
    # Only check columns that exist
    check_cols = [c for c in check_cols if c in merged.columns]

    n_before_drop = len(merged)
    merged = merged.dropna(subset=check_cols).copy()
    n_dropped = n_before_drop - len(merged)
    summary["n_after_missing_drop"] = len(merged)
    log.info(f"  Dropped {n_dropped} samples with missing values → "
             f"{len(merged)} remain")

    if len(merged) == 0:
        log.error(f"  No samples remain after dropping missing for {dataset}")
        summary["status"] = "EMPTY_AFTER_NA_DROP"
        return summary

    # ── 8. Standardize PGS (mean=0, sd=1) within this dataset ───────────
    for pgs_col in available_pgs:
        col_mean = merged[pgs_col].mean()
        col_std = merged[pgs_col].std()
        if col_std > 0:
            merged[pgs_col] = (merged[pgs_col] - col_mean) / col_std
        else:
            log.warning(f"  PGS column {pgs_col} has zero variance in {dataset}")
            merged[pgs_col] = 0.0

    # ── 9. Select output columns and save ────────────────────────────────
    out_cols = (
        ["IID", "phenotype"]
        + [DOSAGE_COLUMN]
        + REQUIRED_COVARIATES
        + available_pgs
    )
    # Add source column if available (useful for downstream meta-analysis)
    if "source" in merged.columns:
        out_cols.append("source")

    out_df = merged[out_cols].copy()

    out_path = outdir / f"{dataset}_{phenotype}_analysis_ready.tsv"
    out_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")

    n_cases = (out_df["phenotype"] == 1).sum()
    n_controls = (out_df["phenotype"] == 0).sum()
    summary["n_cases"] = int(n_cases)
    summary["n_controls"] = int(n_controls)
    log.info(f"  Output: {out_path}")
    log.info(f"  Cases: {n_cases}  Controls: {n_controls}  "
             f"Total: {len(out_df)}")

    # Quick dosage distribution among cases vs controls
    dos_cases = out_df.loc[out_df["phenotype"] == 1, DOSAGE_COLUMN]
    dos_ctrls = out_df.loc[out_df["phenotype"] == 0, DOSAGE_COLUMN]
    log.info(f"  {DOSAGE_COLUMN} — cases:  mean={dos_cases.mean():.4f}  "
             f"sd={dos_cases.std():.4f}")
    log.info(f"  {DOSAGE_COLUMN} — controls: mean={dos_ctrls.mean():.4f}  "
             f"sd={dos_ctrls.std():.4f}")

    return summary


# ── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge covariates, PGS z-scores, and rs55705857 dosage "
                    "into analysis-ready files for ICVF conditioning.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Phenotype definitions:
  idhmt         IDH-mutant cases (idh==1) vs controls (case==0)
  idhmt_intact  IDH-mutant 1p/19q-intact (idh==1, pq==0) vs controls
  idhmt_codel   IDH-mutant 1p/19q-codel  (idh==1, pq==1) vs controls
        """,
    )
    parser.add_argument(
        "--scores-dir", required=True, type=Path,
        help="Directory with {dataset}-covariates.csv and "
             "{dataset}.scores.z-scores.txt.gz files",
    )
    parser.add_argument(
        "--dosage-dir", required=True, type=Path,
        help="Directory with {dataset}_rs55705857_dosage.tsv files "
             "(output of 01_extract_dosage.sh)",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Output directory for analysis-ready TSV files",
    )
    parser.add_argument(
        "--phenotype", default="idhmt",
        choices=sorted(VALID_PHENOTYPES),
        help="Phenotype definition (default: idhmt)",
    )
    parser.add_argument(
        "--datasets", nargs="+", default=DATASETS,
        help=f"Dataset names to process (default: {' '.join(DATASETS)})",
    )
    args = parser.parse_args()

    # Validate inputs
    if not args.scores_dir.is_dir():
        log.error(f"Scores directory not found: {args.scores_dir}")
        sys.exit(1)
    if not args.dosage_dir.is_dir():
        log.error(f"Dosage directory not found: {args.dosage_dir}")
        sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)

    log.info(f"Scores dir:  {args.scores_dir}")
    log.info(f"Dosage dir:  {args.dosage_dir}")
    log.info(f"Output dir:  {args.outdir}")
    log.info(f"Phenotype:   {args.phenotype}")
    log.info(f"Datasets:    {args.datasets}")
    log.info(f"ICVF PGS:    {len(ICVF_PGS_COLUMNS)} models")

    # Process each dataset
    summaries = []
    for dataset in args.datasets:
        summary = process_dataset(
            dataset=dataset,
            scores_dir=args.scores_dir,
            dosage_dir=args.dosage_dir,
            outdir=args.outdir,
            phenotype=args.phenotype,
        )
        summaries.append(summary)

    # ── Write summary table ──────────────────────────────────────────────
    summary_df = pd.DataFrame(summaries)
    summary_path = args.outdir / "merge_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    log.info("=" * 60)
    log.info("MERGE SUMMARY")
    log.info("=" * 60)
    for _, row in summary_df.iterrows():
        log.info(
            f"  {row['dataset']:6s}  "
            f"cases={row['n_cases']:5d}  "
            f"controls={row['n_controls']:5d}  "
            f"total={row['n_after_missing_drop']:5d}  "
            f"status={row['status']}"
        )

    total_cases = summary_df["n_cases"].sum()
    total_controls = summary_df["n_controls"].sum()
    log.info(f"  TOTAL   cases={total_cases:5d}  controls={total_controls:5d}  "
             f"total={total_cases + total_controls:5d}")
    log.info(f"Summary written to: {summary_path}")

    # Exit with error if any dataset failed
    failed = summary_df[summary_df["status"] != "OK"]
    if len(failed) > 0:
        log.warning(f"{len(failed)} dataset(s) had issues: "
                    f"{', '.join(failed['dataset'].tolist())}")
        sys.exit(1)

    log.info("All datasets processed successfully.")


if __name__ == "__main__":
    main()
