#!/usr/bin/env python3
"""
08_tables.py — Generate publication-quality CSV tables.

Tables:
  1. sample_characteristics.csv  — demographics per dataset
  2. primary_results.csv         — wide-format results across all analyses
  3. rs55705857_info.csv         — variant info per dataset
  4. pgs_coverage.csv            — PGS variant coverage and LD-pruning counts

Usage
-----
python 08_tables.py \
    --data-dir merged_data/ \
    --results-dir results/ \
    --outdir tables/
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.analysis_utils import (
    ICVF_PGS_IDS,
    PGS_WITHOUT_RS55705857,
    RS55705857_COL,
    DATASET_NAMES,
    discover_datasets,
    load_dataset,
)

PGS_SHORT_NAMES = {
    "PGS001454": "Body CC",
    "PGS001456": "Cerebral peduncle R",
    "PGS001457": "Cing. cingulate gyrus L",
    "PGS001458": "Cing. cingulate gyrus R",
    "PGS001459": "Cing. hippocampus L",
    "PGS001460": "Cing. hippocampus R",
    "PGS001466": "Genu CC",
    "PGS001471": "Mid. cerebellar peduncle",
    "PGS001474": "Post. int. capsule L",
    "PGS001478": "Retrolent. int. capsule L",
    "PGS001479": "Retrolent. int. capsule R",
    "PGS001480": "Sagittal stratum L",
    "PGS001481": "Sagittal stratum R",
    "PGS001484": "Sup. corona radiata L",
    "PGS001485": "Sup. corona radiata R",
    "PGS001662": "Acoustic radiation R",
    "PGS001669": "Forceps major",
    "PGS001679": "Parahippocampal cing. R",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate publication tables.")
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--phenotype", default="idhmt")
    return parser.parse_args()


def fmt_or_ci(row, or_col="or_", lo_col="or_lower", hi_col="or_upper"):
    o = row.get(or_col, np.nan)
    lo = row.get(lo_col, np.nan)
    hi = row.get(hi_col, np.nan)
    if not np.isfinite(o):
        return "\u2014"
    return f"{o:.2f} ({lo:.2f}\u2013{hi:.2f})"


def fmt_p(val):
    if not np.isfinite(val):
        return "\u2014"
    if val < 0.001:
        return f"{val:.1e}"
    return f"{val:.3f}"


# =========================================================================
# Table 1: Sample characteristics
# =========================================================================

def table1_sample_characteristics(data_dir, phenotype, outdir):
    datasets = discover_datasets(data_dir, phenotype)
    rows = []
    for ds_name, ds_path in sorted(datasets.items()):
        df = load_dataset(ds_path)
        n_total = len(df)
        n_cases = int((df["phenotype"] == 2).sum())
        n_controls = int((df["phenotype"] == 1).sum())

        dosage = df[RS55705857_COL].dropna()
        geno = dosage.round().astype(int)
        n_AA = int((geno == 0).sum())
        n_AG = int((geno == 1).sum())
        n_GG = int((geno == 2).sum())
        carrier_pct = (n_AG + n_GG) / len(geno) * 100 if len(geno) > 0 else np.nan

        age = df["age"].dropna()
        sex = df["sex"].dropna()

        rows.append(dict(
            dataset=ds_name, n_total=n_total, n_cases=n_cases, n_controls=n_controls,
            rs55705857_AA=n_AA, rs55705857_AG=n_AG, rs55705857_GG=n_GG,
            carrier_pct=round(carrier_pct, 1),
            age_mean=round(age.mean(), 1), age_sd=round(age.std(), 1),
            pct_male=round((sex == 1).mean() * 100, 1) if len(sex) > 0 else np.nan,
        ))

    out = pd.DataFrame(rows)
    path = os.path.join(outdir, "sample_characteristics.csv")
    out.to_csv(path, index=False)
    print(f"  Table 1 written: {path} ({len(out)} rows)")


# =========================================================================
# Table 2: Primary results (wide format)
# =========================================================================

def table2_primary_results(results_dir, outdir):
    def _load_meta(fname, model=None, term=None):
        fpath = os.path.join(results_dir, fname)
        if not os.path.isfile(fpath):
            return pd.DataFrame()
        df = pd.read_csv(fpath)
        mask = df["dataset"] == "meta"
        if model:
            mask &= df["model"] == model
        if term:
            mask &= df["term"] == term
        return df.loc[mask].set_index("pgs_id")

    a2_uncond = _load_meta("analysis2_conditional_results.csv", model="unconditional")
    a2_cond = _load_meta("analysis2_conditional_results.csv", model="conditional")
    a1_aa = _load_meta("analysis1_stratified_results.csv", model="stratified_noncarrier")
    a3_inter = _load_meta("analysis3_interaction_results.csv", term="pgs_x_rs55705857")
    a4_pruned = _load_meta("analysis4_ld_pruned_results.csv", model="ld_pruned")

    result_rows = []
    for pgs_id in ICVF_PGS_IDS:
        row = {
            "pgs_id": pgs_id,
            "tract_name": PGS_SHORT_NAMES.get(pgs_id, pgs_id),
            "negative_control": pgs_id in PGS_WITHOUT_RS55705857,
        }

        for label, src in [("uncond", a2_uncond), ("cond", a2_cond),
                           ("strat_AA", a1_aa), ("pruned", a4_pruned)]:
            if pgs_id in src.index:
                r = src.loc[pgs_id]
                row[f"{label}_OR_CI"] = fmt_or_ci(r)
                row[f"{label}_p"] = fmt_p(r.get("p", np.nan))
                row[f"{label}_fdr_q"] = fmt_p(r.get("fdr_q", np.nan))
            else:
                row[f"{label}_OR_CI"] = "\u2014"
                row[f"{label}_p"] = "\u2014"
                row[f"{label}_fdr_q"] = "\u2014"

        if pgs_id in a3_inter.index:
            r = a3_inter.loc[pgs_id]
            b, s = r.get("beta", np.nan), r.get("se", np.nan)
            row["interaction_beta_SE"] = f"{b:.3f} ({s:.3f})" if np.isfinite(b) else "\u2014"
            row["interaction_p"] = fmt_p(r.get("p", np.nan))
            row["interaction_fdr_q"] = fmt_p(r.get("fdr_q", np.nan))
        else:
            row["interaction_beta_SE"] = "\u2014"
            row["interaction_p"] = "\u2014"
            row["interaction_fdr_q"] = "\u2014"

        result_rows.append(row)

    out = pd.DataFrame(result_rows)
    path = os.path.join(outdir, "primary_results.csv")
    out.to_csv(path, index=False)
    print(f"  Table 2 written: {path} ({len(out)} rows)")


# =========================================================================
# Table 3: rs55705857 variant info per dataset
# =========================================================================

def table3_rs55705857_info(data_dir, phenotype, outdir):
    datasets = discover_datasets(data_dir, phenotype)
    rows = []
    for ds_name, ds_path in sorted(datasets.items()):
        df = load_dataset(ds_path)
        dosage = df[RS55705857_COL].dropna()
        af = dosage.mean() / 2 if len(dosage) > 0 else np.nan

        geno = dosage.round().astype(int)
        n_AA = int((geno == 0).sum())
        n_AG = int((geno == 1).sum())
        n_GG = int((geno == 2).sum())
        carrier_pct = (n_AG + n_GG) / len(geno) * 100 if len(geno) > 0 else np.nan

        cases = df[df["phenotype"] == 2]
        controls = df[df["phenotype"] == 1]
        case_carrier = ((cases[RS55705857_COL].round() >= 1).sum() / len(cases) * 100
                        if len(cases) > 0 else np.nan)
        ctrl_carrier = ((controls[RS55705857_COL].round() >= 1).sum() / len(controls) * 100
                        if len(controls) > 0 else np.nan)

        rows.append(dict(
            dataset=ds_name,
            risk_allele_freq=round(af, 4) if np.isfinite(af) else np.nan,
            n_AA=n_AA, n_AG=n_AG, n_GG=n_GG,
            carrier_pct=round(carrier_pct, 1),
            case_carrier_pct=round(case_carrier, 1) if np.isfinite(case_carrier) else np.nan,
            control_carrier_pct=round(ctrl_carrier, 1) if np.isfinite(ctrl_carrier) else np.nan,
        ))

    out = pd.DataFrame(rows)
    path = os.path.join(outdir, "rs55705857_info.csv")
    out.to_csv(path, index=False)
    print(f"  Table 3 written: {path} ({len(out)} rows)")


# =========================================================================
# Table 4: PGS coverage and LD pruning
# =========================================================================

def table4_pgs_coverage(results_dir, outdir):
    pruning_summary = None
    for candidate in [
        os.path.join(results_dir, "..", "ld_pruned_scores", "pruning_summary.tsv"),
        os.path.join(results_dir, "pruning_summary.tsv"),
    ]:
        if os.path.isfile(candidate):
            pruning_summary = pd.read_csv(candidate, sep="\t")
            break

    rows = []
    for pgs_id in ICVF_PGS_IDS:
        row = {
            "pgs_id": pgs_id,
            "tract_name": PGS_SHORT_NAMES.get(pgs_id, pgs_id),
            "rs55705857_in_scoring": pgs_id not in PGS_WITHOUT_RS55705857,
        }
        if pruning_summary is not None and pgs_id in pruning_summary["pgs_id"].values:
            ps = pruning_summary[pruning_summary["pgs_id"] == pgs_id].iloc[0]
            row["n_total_variants"] = int(ps.get("n_total", 0))
            row["n_pruned"] = int(ps.get("n_pruned", 0))
            row["n_remaining"] = int(ps.get("n_remaining", 0))
            row["pct_removed"] = round(ps.get("pct_removed", 0), 2)
        else:
            row["n_total_variants"] = np.nan
            row["n_pruned"] = np.nan
            row["n_remaining"] = np.nan
            row["pct_removed"] = np.nan
        rows.append(row)

    out = pd.DataFrame(rows)
    path = os.path.join(outdir, "pgs_coverage.csv")
    out.to_csv(path, index=False)
    print(f"  Table 4 written: {path} ({len(out)} rows)")


# =========================================================================
# Main
# =========================================================================

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("=== Generating publication tables ===")

    print("[Table 1] Sample characteristics")
    table1_sample_characteristics(args.data_dir, args.phenotype, args.outdir)

    print("[Table 2] Primary results (wide format)")
    table2_primary_results(args.results_dir, args.outdir)

    print("[Table 3] rs55705857 variant info")
    table3_rs55705857_info(args.data_dir, args.phenotype, args.outdir)

    print("[Table 4] PGS coverage and LD pruning")
    table4_pgs_coverage(args.results_dir, args.outdir)

    print("\n=== All tables complete ===")


if __name__ == "__main__":
    main()
