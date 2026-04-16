#!/usr/bin/env python3
"""
params.py — Shared configuration for the ICVF conditioning analysis pipeline.

Centralises every magic constant (PGS IDs, variant coordinates, covariate
lists, phenotype definitions) so that the individual analysis scripts import
from one place and stay in sync.

Usage
-----
    from utils.params import (
        ICVF_PGS_IDS, PGS_SHORT_NAMES, NEGATIVE_CONTROL_PGS,
        RS55705857_HG19_POS, RS55705857_ALLELES,
        STANDARD_COVARIATES, DATASETS_DEFAULT, PHENOTYPE_DEFINITIONS,
    )

Or load everything as a dict:

    from utils.params import load_params
    cfg = load_params()                     # all defaults
    cfg = load_params("my_config.tsv")      # override with a file
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

# ═══════════════════════════════════════════════════════════════════════
#  rs55705857 — the major IDH-mutant glioma risk variant (OR ≈ 6.5)
# ═══════════════════════════════════════════════════════════════════════

RS55705857_HG19_POS = "8:130645692"       # no chr prefix, matches VCF convention
RS55705857_HG38_POS = "chr8:129633446"    # UCSC-style with chr prefix
RS55705857_ALLELES = {
    "ref": "A",
    "alt": "G",
    "risk": "G",
}

# ═══════════════════════════════════════════════════════════════════════
#  ICVF polygenic scores — 18 white-matter tract PGS from PGS Catalog
# ═══════════════════════════════════════════════════════════════════════

ICVF_PGS_IDS: list[str] = [
    "PGS001454",   # Body of corpus callosum
    "PGS001456",   # Cerebral peduncle (R)
    "PGS001457",   # Cingulum — cingulate gyrus (L)
    "PGS001458",   # Cingulum — cingulate gyrus (R)
    "PGS001459",   # Cingulum — hippocampus (L)
    "PGS001460",   # Cingulum — hippocampus (R)
    "PGS001466",   # Genu of corpus callosum
    "PGS001471",   # Middle cerebellar peduncle
    "PGS001474",   # Posterior limb of internal capsule (L)
    "PGS001478",   # Retrolenticular part of internal capsule (L)
    "PGS001479",   # Retrolenticular part of internal capsule (R)
    "PGS001480",   # Sagittal stratum (L)
    "PGS001481",   # Sagittal stratum (R)
    "PGS001484",   # Superior corona radiata (L)
    "PGS001485",   # Superior corona radiata (R)
    "PGS001662",   # Acoustic radiation (R)          [WA-ICVF]
    "PGS001669",   # Forceps major                   [WA-ICVF]
    "PGS001679",   # Parahippocampal cingulum (R)    [WA-ICVF]
]

# Full official trait names from the PGS Catalog
PGS_FULL_NAMES: dict[str, str] = {
    "PGS001454": "Mean ICVF in body of corpus callosum on FA skeleton",
    "PGS001456": "Mean ICVF in cerebral peduncle on FA skeleton (R)",
    "PGS001457": "Mean ICVF in cingulum cingulate gyrus on FA skeleton (L)",
    "PGS001458": "Mean ICVF in cingulum cingulate gyrus on FA skeleton (R)",
    "PGS001459": "Mean ICVF in cingulum hippocampus on FA skeleton (L)",
    "PGS001460": "Mean ICVF in cingulum hippocampus on FA skeleton (R)",
    "PGS001466": "Mean ICVF in genu of corpus callosum on FA skeleton",
    "PGS001471": "Mean ICVF in middle cerebellar peduncle on FA skeleton",
    "PGS001474": "Mean ICVF in posterior limb of internal capsule on FA skeleton (L)",
    "PGS001478": "Mean ICVF in retrolenticular part of internal capsule on FA skeleton (L)",
    "PGS001479": "Mean ICVF in retrolenticular part of internal capsule on FA skeleton (R)",
    "PGS001480": "Mean ICVF in sagittal stratum on FA skeleton (L)",
    "PGS001481": "Mean ICVF in sagittal stratum on FA skeleton (R)",
    "PGS001484": "Mean ICVF in superior corona radiata on FA skeleton (L)",
    "PGS001485": "Mean ICVF in superior corona radiata on FA skeleton (R)",
    "PGS001662": "WA ICVF in tract acoustic radiation (R)",
    "PGS001669": "WA ICVF in tract forceps major",
    "PGS001679": "WA ICVF in tract parahippocampal part of cingulum (R)",
}

# Short names for figures and table labels
PGS_SHORT_NAMES: dict[str, str] = {
    "PGS001454": "Body CC",
    "PGS001456": "Cerebral peduncle (R)",
    "PGS001457": "Cingulum CG (L)",
    "PGS001458": "Cingulum CG (R)",
    "PGS001459": "Cingulum Hipp (L)",
    "PGS001460": "Cingulum Hipp (R)",
    "PGS001466": "Genu CC",
    "PGS001471": "Middle cerebellar ped",
    "PGS001474": "PLIC (L)",
    "PGS001478": "RLIC (L)",
    "PGS001479": "RLIC (R)",
    "PGS001480": "Sagittal stratum (L)",
    "PGS001481": "Sagittal stratum (R)",
    "PGS001484": "Sup corona radiata (L)",
    "PGS001485": "Sup corona radiata (R)",
    "PGS001662": "Acoustic radiation (R)",
    "PGS001669": "Forceps major",
    "PGS001679": "Parahipp cingulum (R)",
}

# ── Negative control PGS ─────────────────────────────────────────────
# These 3 ICVF PGS do NOT contain rs55705857 in their scoring file,
# so they serve as negative controls: conditioning on rs55705857 should
# NOT change their association with IDHmt risk.
NEGATIVE_CONTROL_PGS: list[str] = [
    "PGS001456",   # Cerebral peduncle (R)
    "PGS001662",   # Acoustic radiation (R)
    "PGS001679",   # Parahippocampal cingulum (R)
]

# ═══════════════════════════════════════════════════════════════════════
#  Covariates and dataset defaults
# ═══════════════════════════════════════════════════════════════════════

STANDARD_COVARIATES: list[str] = [
    "age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
]

DATASETS_DEFAULT: list[str] = ["cidr", "i370", "onco", "tcga"]

# ═══════════════════════════════════════════════════════════════════════
#  Phenotype definitions
# ═══════════════════════════════════════════════════════════════════════
#
# Each entry maps phenotype name → dict with:
#   label       : human-readable label
#   case_filter : dict of column→value(s) for case selection
#   ctrl_filter : dict of column→value(s) for control selection
#   description : brief description for logs / reports
#
# The merge script applies these filters; analysis scripts use the
# resulting binary 'phenotype' column (1 = case, 0 = control).

PHENOTYPE_DEFINITIONS: dict[str, dict[str, Any]] = {
    "idhmt": {
        "label": "IDH-mutant glioma",
        "case_filter": {"idh": [1]},
        "ctrl_filter": {"case": [0]},
        "description": "IDH-mutant glioma cases vs population controls",
    },
    "idhmt_intact": {
        "label": "IDH-mutant, 1p/19q-intact",
        "case_filter": {"idh": [1], "pq": [0]},
        "ctrl_filter": {"case": [0]},
        "description": (
            "IDH-mutant, 1p/19q-intact (astrocytoma) cases vs population controls"
        ),
    },
    "idhmt_codel": {
        "label": "IDH-mutant, 1p/19q-codeleted",
        "case_filter": {"idh": [1], "pq": [1]},
        "ctrl_filter": {"case": [0]},
        "description": (
            "IDH-mutant, 1p/19q-codeleted (oligodendroglioma) cases vs "
            "population controls"
        ),
    },
}

# ═══════════════════════════════════════════════════════════════════════
#  Configuration loader
# ═══════════════════════════════════════════════════════════════════════


def load_params(config_file: str | Path | None = None) -> dict[str, Any]:
    """Return pipeline parameters as a single dict.

    If *config_file* is ``None``, returns hard-coded defaults from this
    module.  If a path is given, reads a two-column TSV (key<TAB>value)
    and overlays those values on top of the defaults.  This allows the
    ``run_pipeline.sh`` master script to serialise runtime parameters
    into a file that individual Python steps can read.

    Returns
    -------
    dict with keys:
        icvf_pgs_ids, pgs_short_names, pgs_full_names,
        negative_control_pgs, rs55705857_hg19_pos, rs55705857_hg38_pos,
        rs55705857_alleles, standard_covariates, datasets,
        phenotype, phenotype_definitions, output_dir, ...
        (plus any extra keys from config file)
    """
    params: dict[str, Any] = {
        # PGS metadata
        "icvf_pgs_ids": ICVF_PGS_IDS,
        "pgs_short_names": PGS_SHORT_NAMES,
        "pgs_full_names": PGS_FULL_NAMES,
        "negative_control_pgs": NEGATIVE_CONTROL_PGS,
        # Variant
        "rs55705857_hg19_pos": RS55705857_HG19_POS,
        "rs55705857_hg38_pos": RS55705857_HG38_POS,
        "rs55705857_alleles": RS55705857_ALLELES,
        # Covariates and datasets
        "standard_covariates": STANDARD_COVARIATES,
        "datasets": DATASETS_DEFAULT,
        # Phenotype
        "phenotype": "idhmt",
        "phenotype_definitions": PHENOTYPE_DEFINITIONS,
        # Paths (overridable from config file)
        "output_dir": "./output",
    }

    if config_file is not None:
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        with open(config_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("key\t"):
                    continue
                parts = line.split("\t", 1)
                key = parts[0]
                val = parts[1] if len(parts) > 1 else ""
                # Convert comma-separated lists back to lists
                if key == "datasets":
                    params[key] = [v.strip() for v in val.split(",") if v.strip()]
                else:
                    params[key] = val

    return params
