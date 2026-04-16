#!/usr/bin/env python3
"""
Shared utilities for ICVF conditioning analysis pipeline.

Provides logistic regression wrappers, inverse-variance weighted
fixed-effects meta-analysis, FDR correction, and data loading helpers.
"""

import sys
import logging
import warnings
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STANDARD_COVARIATES = [
    "age", "sex",
    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
]

ICVF_PGS_IDS = [
    "PGS001454", "PGS001456", "PGS001457", "PGS001458", "PGS001459",
    "PGS001460", "PGS001466", "PGS001471", "PGS001474", "PGS001478",
    "PGS001479", "PGS001480", "PGS001481", "PGS001484", "PGS001485",
    "PGS001662", "PGS001669", "PGS001679",
]

# PGS IDs whose scoring files do NOT contain rs55705857
PGS_WITHOUT_RS55705857 = {"PGS001456", "PGS001662", "PGS001679"}

RS55705857_COL = "rs55705857_dosage"

# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------

def setup_logging(name: str = "icvf_analysis", level: int = logging.INFO):
    """Configure structured stderr logging."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        fmt = logging.Formatter(
            "[%(asctime)s] %(levelname)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dataset(filepath: str) -> pd.DataFrame:
    """Load an analysis-ready TSV file produced by 02_merge_data.py.

    Parameters
    ----------
    filepath : str
        Path to the TSV file.

    Returns
    -------
    pd.DataFrame
    """
    df = pd.read_csv(filepath, sep="\t")
    return df

# ---------------------------------------------------------------------------
# Logistic regression
# ---------------------------------------------------------------------------

def run_logistic(
    df: pd.DataFrame,
    pgs_col: str,
    covariates: list[str],
    phenotype_col: str = "phenotype",
) -> dict:
    """Run logistic regression: phenotype ~ pgs_col + covariates.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain *phenotype_col*, *pgs_col*, and all *covariates*.
    pgs_col : str
        Name of the ICVF PGS column (predictor of interest).
    covariates : list[str]
        Additional covariate column names (e.g. age, sex, PCs).
    phenotype_col : str
        Column with case/control coding (2 = case, 1 = control).

    Returns
    -------
    dict
        Keys: beta, se, z, p, n, n_cases, n_controls, or_, or_lower,
        or_upper, converged.  Values are NaN when convergence fails.
    """
    nan_result = dict(
        beta=np.nan, se=np.nan, z=np.nan, p=np.nan,
        n=np.nan, n_cases=np.nan, n_controls=np.nan,
        or_=np.nan, or_lower=np.nan, or_upper=np.nan,
        converged=False,
    )

    # Subset to required columns and drop missing
    cols = [phenotype_col, pgs_col] + list(covariates)
    sub = df[cols].dropna()

    # Recode phenotype to 1/0. Handle both conventions:
    #   PLINK: 2=case, 1=control  →  map {2:1, 1:0}
    #   Standard: 1=case, 0=control  →  already correct
    raw_vals = set(sub[phenotype_col].dropna().unique())
    if raw_vals <= {2, 1}:
        y = sub[phenotype_col].map({2: 1, 1: 0})
    elif raw_vals <= {1, 0}:
        y = sub[phenotype_col].astype(float)
    else:
        raise ValueError(f"Unexpected phenotype values: {raw_vals}")
    if y.isna().any():
        # Unexpected phenotype values
        sub = sub.loc[y.notna()]
        y = y.loc[y.notna()]

    n_cases = int(y.sum())
    n_controls = int(len(y) - n_cases)
    n = len(y)

    if n_cases < 5 or n_controls < 5:
        nan_result.update(n=n, n_cases=n_cases, n_controls=n_controls)
        return nan_result

    # Design matrix: PGS column first, then covariates, plus intercept
    X = sm.add_constant(sub[[pgs_col] + list(covariates)].astype(float))

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = Logit(y.astype(float), X)
            result = model.fit(disp=0, maxiter=100, method="newton", warn_convergence=False)

            # Check convergence via MLE return code
            converged = result.mle_retvals.get("converged", False) if hasattr(result, "mle_retvals") else True

            # Extract the coefficient for the PGS column (index 1, after const)
            idx = 1  # PGS is the first predictor after the constant
            beta = result.params.iloc[idx]
            se = result.bse.iloc[idx]
            z_val = result.tvalues.iloc[idx]
            p_val = result.pvalues.iloc[idx]

            or_val = np.exp(beta)
            or_lower = np.exp(beta - 1.96 * se)
            or_upper = np.exp(beta + 1.96 * se)

            return dict(
                beta=beta, se=se, z=z_val, p=p_val,
                n=n, n_cases=n_cases, n_controls=n_controls,
                or_=or_val, or_lower=or_lower, or_upper=or_upper,
                converged=converged,
            )

    except Exception:
        nan_result.update(n=n, n_cases=n_cases, n_controls=n_controls)
        return nan_result


def run_logistic_interaction(
    df: pd.DataFrame,
    pgs_col: str,
    interaction_col: str,
    covariates: list[str],
    phenotype_col: str = "phenotype",
) -> dict:
    """Logistic regression with an interaction term.

    Model: phenotype ~ pgs_col + interaction_col + pgs_col*interaction_col + covariates

    Returns a dict with results for *main* (pgs_col), *interaction_var*
    (interaction_col), and *interaction_term* (pgs_col × interaction_col).
    """
    nan_single = dict(beta=np.nan, se=np.nan, z=np.nan, p=np.nan)
    nan_result = dict(
        main=nan_single.copy(),
        interaction_var=nan_single.copy(),
        interaction_term=nan_single.copy(),
        n=np.nan, n_cases=np.nan, n_controls=np.nan, converged=False,
    )

    cols = [phenotype_col, pgs_col, interaction_col] + list(covariates)
    sub = df[cols].dropna()

    # Recode phenotype (auto-detect 2/1 vs 1/0 convention)
    raw_vals = set(sub[phenotype_col].dropna().unique())
    if raw_vals <= {2, 1}:
        y = sub[phenotype_col].map({2: 1, 1: 0})
    elif raw_vals <= {1, 0}:
        y = sub[phenotype_col].astype(float)
    else:
        raise ValueError(f"Unexpected phenotype values: {raw_vals}")
    if y.isna().any():
        sub = sub.loc[y.notna()]
        y = y.loc[y.notna()]

    n_cases = int(y.sum())
    n_controls = int(len(y) - n_cases)
    n = len(y)

    if n_cases < 5 or n_controls < 5:
        nan_result.update(n=n, n_cases=n_cases, n_controls=n_controls)
        return nan_result

    # Build design matrix with interaction term
    interaction_name = f"{pgs_col}_x_{interaction_col}"
    sub = sub.copy()
    sub[interaction_name] = sub[pgs_col].astype(float) * sub[interaction_col].astype(float)

    predictor_cols = [pgs_col, interaction_col, interaction_name] + list(covariates)
    X = sm.add_constant(sub[predictor_cols].astype(float))

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = Logit(y.astype(float), X)
            result = model.fit(disp=0, maxiter=100, method="newton", warn_convergence=False)
            converged = result.mle_retvals.get("converged", False) if hasattr(result, "mle_retvals") else True

            def _extract(col_name):
                idx = list(X.columns).index(col_name)
                b = result.params.iloc[idx]
                s = result.bse.iloc[idx]
                return dict(
                    beta=b, se=s,
                    z=result.tvalues.iloc[idx],
                    p=result.pvalues.iloc[idx],
                    or_=np.exp(b),
                    or_lower=np.exp(b - 1.96 * s),
                    or_upper=np.exp(b + 1.96 * s),
                )

            return dict(
                main=_extract(pgs_col),
                interaction_var=_extract(interaction_col),
                interaction_term=_extract(interaction_name),
                n=n, n_cases=n_cases, n_controls=n_controls,
                converged=converged,
            )

    except Exception:
        nan_result.update(n=n, n_cases=n_cases, n_controls=n_controls)
        return nan_result


# ---------------------------------------------------------------------------
# Meta-analysis
# ---------------------------------------------------------------------------

def meta_analyze_ivw(results_list: list[dict]) -> dict:
    """Inverse-variance weighted fixed-effects meta-analysis.

    Parameters
    ----------
    results_list : list[dict]
        Each dict must have keys *beta* and *se*.  Entries with NaN beta/se
        are silently skipped.

    Returns
    -------
    dict
        meta_beta, meta_se, meta_z, meta_p, meta_or, meta_or_lower,
        meta_or_upper, Q, Q_p, I2, n_datasets.
    """
    nan_result = dict(
        meta_beta=np.nan, meta_se=np.nan, meta_z=np.nan, meta_p=np.nan,
        meta_or=np.nan, meta_or_lower=np.nan, meta_or_upper=np.nan,
        Q=np.nan, Q_p=np.nan, I2=np.nan, n_datasets=0,
    )

    # Filter valid entries
    valid = [
        r for r in results_list
        if np.isfinite(r.get("beta", np.nan)) and np.isfinite(r.get("se", np.nan)) and r["se"] > 0
    ]

    n_datasets = len(valid)
    if n_datasets == 0:
        return nan_result

    betas = np.array([r["beta"] for r in valid])
    ses = np.array([r["se"] for r in valid])
    weights = 1.0 / (ses ** 2)

    meta_beta = np.sum(weights * betas) / np.sum(weights)
    meta_se = np.sqrt(1.0 / np.sum(weights))
    meta_z = meta_beta / meta_se
    meta_p = 2 * sp_stats.norm.sf(np.abs(meta_z))

    meta_or = np.exp(meta_beta)
    meta_or_lower = np.exp(meta_beta - 1.96 * meta_se)
    meta_or_upper = np.exp(meta_beta + 1.96 * meta_se)

    # Cochran's Q and I²
    Q = np.sum(weights * (betas - meta_beta) ** 2)
    Q_df = max(n_datasets - 1, 1)
    Q_p = 1 - sp_stats.chi2.cdf(Q, Q_df) if n_datasets > 1 else np.nan
    I2 = max(0, (Q - Q_df) / Q * 100) if Q > 0 and n_datasets > 1 else 0.0

    # Aggregate sample counts across datasets
    total_n = sum(r.get("n", 0) for r in results_list if r.get("n") is not None)
    total_cases = sum(r.get("n_cases", 0) for r in results_list if r.get("n_cases") is not None)
    total_controls = sum(r.get("n_controls", 0) for r in results_list if r.get("n_controls") is not None)

    return dict(
        meta_beta=meta_beta, meta_se=meta_se, meta_z=meta_z, meta_p=meta_p,
        meta_or=meta_or, meta_or_lower=meta_or_lower, meta_or_upper=meta_or_upper,
        Q=Q, Q_p=Q_p, I2=I2, n_datasets=n_datasets,
        n=total_n, n_cases=total_cases, n_controls=total_controls,
    )


# ---------------------------------------------------------------------------
# FDR correction
# ---------------------------------------------------------------------------

def apply_fdr(p_values) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like
        Raw p-values.  NaN entries are preserved as NaN in the output.

    Returns
    -------
    np.ndarray
        FDR-adjusted q-values.
    """
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return q

    p_valid = p[mask]
    n = len(p_valid)
    sorted_idx = np.argsort(p_valid)
    sorted_p = p_valid[sorted_idx]

    # BH procedure
    adjusted = sorted_p * n / (np.arange(1, n + 1))
    # Enforce monotonicity (from largest to smallest)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    q_valid = np.empty(n)
    q_valid[sorted_idx] = adjusted
    q[mask] = q_valid
    return q


# ---------------------------------------------------------------------------
# Dataset discovery
# ---------------------------------------------------------------------------

DATASET_NAMES = ["cidr", "i370", "onco", "tcga"]


def discover_datasets(data_dir: str, phenotype: str = "idhmt") -> dict[str, str]:
    """Find per-dataset analysis-ready TSV files in *data_dir*.

    Expected naming: {dataset}_{phenotype}_analysis_ready.tsv

    Returns
    -------
    dict
        Mapping dataset name → file path for files that exist.
    """
    found = {}
    for ds in DATASET_NAMES:
        fname = f"{ds}_{phenotype}_analysis_ready.tsv"
        fpath = os.path.join(data_dir, fname)
        if os.path.isfile(fpath):
            found[ds] = fpath
    return found


# ---------------------------------------------------------------------------
# Result formatting helpers
# ---------------------------------------------------------------------------

def format_result_row(
    pgs_id: str,
    dataset: str,
    model: str,
    res: dict,
    term: str = "pgs",
) -> dict:
    """Build a flat dict suitable for a results DataFrame row."""
    return dict(
        pgs_id=pgs_id,
        dataset=dataset,
        model=model,
        term=term,
        beta=res.get("beta", np.nan),
        se=res.get("se", np.nan),
        z=res.get("z", np.nan),
        p=res.get("p", np.nan),
        or_=res.get("or_", np.nan),
        or_lower=res.get("or_lower", np.nan),
        or_upper=res.get("or_upper", np.nan),
        n=res.get("n", np.nan),
        n_cases=res.get("n_cases", np.nan),
        n_controls=res.get("n_controls", np.nan),
        converged=res.get("converged", np.nan),
    )


def format_meta_row(
    pgs_id: str,
    model: str,
    meta: dict,
    term: str = "pgs",
) -> dict:
    """Build a flat dict for a meta-analysis results row."""
    return dict(
        pgs_id=pgs_id,
        dataset="meta",
        model=model,
        term=term,
        beta=meta.get("meta_beta", np.nan),
        se=meta.get("meta_se", np.nan),
        z=meta.get("meta_z", np.nan),
        p=meta.get("meta_p", np.nan),
        or_=meta.get("meta_or", np.nan),
        or_lower=meta.get("meta_or_lower", np.nan),
        or_upper=meta.get("meta_or_upper", np.nan),
        n=np.nan,
        n_cases=np.nan,
        n_controls=np.nan,
        converged=np.nan,
        Q=meta.get("Q", np.nan),
        Q_p=meta.get("Q_p", np.nan),
        I2=meta.get("I2", np.nan),
        n_datasets=meta.get("n_datasets", 0),
    )


import os  # noqa: E402  (imported at top for type hints but re-stated for clarity)
