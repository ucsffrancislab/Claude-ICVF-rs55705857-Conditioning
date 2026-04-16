# Claude-ICVF-rs55705857-Conditioning


It's NOT an extension of any single prior pipeline:

Prior Pipeline	Relationship
GWAS meta-analysis pipeline (session 7ce561e0)	Produced the subtype-specific summary statistics used in MR. This new pipeline doesn't do GWAS — it uses individual-level data directly.
MR pipeline (session d8e04e5f / your repo)	Tests causal direction (ICVF → glioma) using summary statistics. Completely different statistical framework (2SLS vs. logistic regression).
Multi-omic integration (sessions bb74746a, 70cae33c)	Produced the master cohort, PGS associations, HLA analyses. Closest relative, but that was a broad survey — this is a focused hypothesis test.
HLA risk/survival pipeline (session 9702f1a1)	Different biology entirely.
What it depends on (inputs from prior work):

Individual-level PGS scores (from your UMich scoring)
Covariate files (from the datasets used across all analyses)
rs55705857 dosage (from imputed VCFs)
The finding that 18 ICVF PGS were significant for IDHwt but not IDHmt (from the PGS catalog-wide analysis, which motivated the hypothesis)
The MR results provide complementary evidence but aren't a data dependency
Why it stands alone:

It has a single, specific hypothesis being tested
It uses a distinct statistical approach (conditional/stratified logistic regression + interaction modeling)
It produces its own self-contained results and interpretation
Someone could understand and reproduce it without needing to re-run any prior pipeline — they just need the input data files




Check that we have everything needed.

What it checks:
Check	What it does
1. ICVF PGS in scores	Searches your *.scores.z-scores.txt.gz files for all 18 PGS IDs as column names, counts samples per file
2. rs55705857 in VCFs	Queries chr8:128748020 (hg38) and chr8:130194364 (hg19) via tabix or grep, reports REF/ALT, imputation R², and mean dosage for allele frequency sanity check
2b. Raw genotypes	If provided, checks .bim/.pvar files for the SNP — tells us if it was directly genotyped (best) vs. imputed (still fine)
3. Variant match rates	Searches scoring logs for how many of each PGS model's variants actually matched your data. If no logs, does a quick sanity check on score distributions
4. Covariates	Confirms each dataset has IDH, 1p/19q, age, sex, PCs


```bash

module load htslib

# Minimal (required args only):
bash validate_icvf_conditioning_inputs.sh \
  --scores-dir /path/to/your/scores \
  --vcf-dir /path/to/hg38/vcfs

# Full (all optional checks too):
bash validate_icvf_conditioning_inputs.sh \
  --scores-dir /path/to/your/scores \
  --vcf-dir /path/to/hg38/vcfs \
  --vcf-hg19-dir /path/to/hg19/vcfs \
  --raw-geno-dir /path/to/raw/plink/files \
  --scoring-log-dir /path/to/plink2/scoring/logs \
  --output validation_report.txt

```



# ICVF Conditioning Analysis Pipeline

**Testing whether brain white-matter polygenic scores associate with IDH-mutant glioma risk after conditioning on the major risk variant rs55705857.**

rs55705857 (chr8:130645692, hg19; A>G, risk allele G) confers ~6.5-fold increased risk for IDH-mutant glioma and lies near *MYC*/*CCDC26*. Many ICVF (intra-cellular volume fraction) polygenic scores include this variant or variants in LD with it. This pipeline tests whether the ICVF–glioma association is entirely driven by rs55705857 or whether an independent polygenic signal persists after conditioning.

## Citation

If you use this pipeline, please cite:

> [Your manuscript reference here]

## Requirements

### System tools

| Tool | Version | Required for |
|------|---------|-------------|
| `bash` | ≥ 4.0 | All steps |
| `bcftools` | ≥ 1.15 | Step 1 (dosage extraction); falls back to `tabix + awk` |
| `tabix` | ≥ 1.15 | Step 1 (fallback if bcftools unavailable) |
| `plink2` | ≥ 2.0 | Step 6 only (LD pruning); optional |
| `pdflatex` | any | Step 9 only (report compilation); optional |

### Python (≥ 3.9)

```
pandas >= 1.4
numpy >= 1.21
statsmodels >= 0.13
scipy >= 1.7
matplotlib >= 3.5
seaborn >= 0.12
```

Install with:
```bash
pip install pandas numpy statsmodels scipy matplotlib seaborn
```

## Input Data

### 1. Imputed genotype VCFs (`--vcf-hg19-dir`)

Directory with per-dataset subdirectories, each containing chromosome 8 imputed dosages:

```
vcf_hg19_dir/
├── cidr/
│   ├── chr8.dose.vcf.gz
│   └── chr8.dose.vcf.gz.tbi
├── i370/
│   ├── chr8.dose.vcf.gz
│   └── chr8.dose.vcf.gz.tbi
├── onco/
│   └── ...
└── tcga/
    └── ...
```

**Key format requirements:**
- Chromosome encoded as `8` (no `chr` prefix)
- VCF ID column uses `CHROM:POS` format (e.g., `8:130645692`)
- `DS` (dosage) field present in FORMAT column
- Gzipped with tabix `.tbi` index

### 2. PGS scores and covariates (`--scores-dir`)

```
scores_dir/
├── cidr-covariates.csv       # Covariate CSV per dataset
├── cidr.scores.z-scores.txt.gz  # PGS z-scores (tab-sep, gzipped)
├── i370-covariates.csv
├── i370.scores.z-scores.txt.gz
├── onco-covariates.csv
├── onco.scores.z-scores.txt.gz
├── tcga-covariates.csv
└── tcga.scores.z-scores.txt.gz
```

**Covariate CSV columns:**

| Column | Type | Description |
|--------|------|-------------|
| `IID` | string | Sample identifier |
| `case` | 0/1 | 0 = population control, 1 = glioma case |
| `idh` | 0/1/9/NaN | 0 = wildtype, 1 = mutant, 9 = unknown |
| `pq` | 0/1/NaN | 0 = 1p/19q intact, 1 = co-deleted |
| `age` | numeric | Age at diagnosis/enrollment |
| `sex` | M/F | Biological sex |
| `PC1`–`PC8` | numeric | Ancestry principal components |
| `exclude` | 0/1 | 1 = exclude from analysis |
| `source` | string | Dataset source label |

**PGS z-score files:**
- Tab-separated, gzipped
- First column: sample ID (`sample` or `IID`)
- Contains columns for all scored PGS models (pipeline uses 18 ICVF-related columns)
- Scores should be computed using `pgs-calc` on hg19 imputed VCFs filtered at R² ≥ 0.3

### 3. LD reference panel (`--ld-ref`, optional)

Plink2-format files (`.pgen`/`.pvar`/`.psam` or `.bed`/`.bim`/`.fam`) of a reference population (e.g., 1000 Genomes EUR) for LD-based pruning in Analysis 4.

### 4. PGS Catalog scoring files (`--pgs-catalog-dir`, optional)

Original PGS Catalog scoring files (`.txt.gz`) needed for Analysis 4 to identify and remove variants in LD with rs55705857. Required only when `--ld-ref` is provided.

## Quick Start

Minimal run (Analyses 1–3 only, no LD pruning):

```bash
bash run_pipeline.sh \
    --scores-dir     /data/scores \
    --vcf-hg19-dir   /data/imputed_hg19 \
    --output-dir     ./output \
    --phenotype      idhmt
```

Full run including LD-pruned analysis:

```bash
bash run_pipeline.sh \
    --scores-dir      /data/scores \
    --vcf-hg19-dir    /data/imputed_hg19 \
    --ld-ref          /ref/1000G_EUR/all_phase3 \
    --pgs-catalog-dir /ref/pgs_catalog_files \
    --output-dir      ./output \
    --phenotype       idhmt
```

## Full Usage

```
bash run_pipeline.sh [OPTIONS]

Required:
  --scores-dir DIR        Covariates + PGS z-scores
  --vcf-hg19-dir DIR      Imputed VCFs (or use --skip-extraction)

Optional:
  --ld-ref PATH           LD reference for Analysis 4
  --pgs-catalog-dir DIR   PGS Catalog scoring files (with --ld-ref)
  --output-dir DIR        Output root (default: ./output)
  --phenotype STR         idhmt | idhmt_intact | idhmt_codel
  --datasets LIST         Comma-separated (default: cidr,i370,onco,tcga)
  --skip-extraction       Skip dosage extraction
  --skip-analyses         Skip to figures/tables only

SLURM:
  --slurm                 Submit as batch job
  --slurm-partition STR   Partition (default: normal)
  --slurm-time STR        Wall time (default: 12:00:00)
  --slurm-mem STR         Memory (default: 32G)
  --slurm-cpus INT        CPUs (default: 4)
```

## Output Description

```
output/{phenotype}/
├── dosage/
│   ├── {dataset}_rs55705857_dosage.tsv   # Per-sample dosage
│   └── extraction_summary.tsv            # QC: AF, R², alleles per dataset
├── merged/
│   ├── {dataset}_analysis_ready.tsv      # Merged analysis file
│   └── merge_summary.tsv                 # Sample counts per dataset
├── results/
│   ├── analysis1_stratified_results.csv  # PGS associations stratified by carrier status
│   ├── analysis2_conditional_results.csv # PGS associations conditioned on rs55705857
│   ├── analysis3_interaction_results.csv # PGS × rs55705857 interaction tests
│   ├── analysis4_ld_pruned_results.csv   # LD-pruned PGS associations (if run)
│   └── ld_pruned/                        # LD-pruned scoring files
├── figures/                              # Publication-ready figures
├── tables/                               # Formatted results tables
├── logs/                                 # Per-step log files
└── pipeline_config.tsv                   # Runtime configuration
```

### Analysis-ready files

Each `{dataset}_analysis_ready.tsv` contains:

| Column | Description |
|--------|-------------|
| `IID` | Sample identifier |
| `phenotype` | 1 = case, 0 = control |
| `rs55705857_dosage` | ALT allele (G) dosage [0, 2] |
| `age` | Age |
| `sex` | 0 = female, 1 = male |
| `PC1`–`PC8` | Ancestry PCs |
| `PGS001454`–`PGS001679` | 18 ICVF PGS (standardised: mean=0, sd=1 within dataset) |

## Analysis Overview

### Analysis 1: Stratified by carrier status

Splits samples into rs55705857 carriers (dosage ≥ 0.5) and non-carriers (dosage < 0.5). Tests each ICVF PGS for association with IDHmt risk within each stratum using logistic regression (PGS + covariates). **Tests whether the PGS effect is present in non-carriers** (where rs55705857 does not contribute).

### Analysis 2: Conditional on rs55705857

Fits logistic regression with both PGS and rs55705857 dosage as covariates:

```
phenotype ~ PGS + rs55705857_dosage + age + sex + PC1-PC8
```

Compares to a base model without PGS. **Tests whether PGS adds predictive value beyond rs55705857.** This is the primary analysis.

### Analysis 3: PGS × rs55705857 interaction

Extends the conditional model with an interaction term:

```
phenotype ~ PGS + rs55705857_dosage + PGS × rs55705857_dosage + covariates
```

**Tests whether the PGS effect differs between carriers and non-carriers** (effect modification).

### Analysis 4: LD-pruned PGS (optional)

Removes all variants within an LD window of rs55705857 (r² > 0.1 in the reference panel) from each PGS scoring file, re-scores samples, and repeats the association analysis. **Provides the most stringent test** by physically removing the variant and its proxies from the score.

## 18 ICVF Polygenic Scores

| PGS ID | White-matter tract |
|--------|-------------------|
| PGS001454 | Body of corpus callosum |
| PGS001456 | Cerebral peduncle (R) † |
| PGS001457 | Cingulum — cingulate gyrus (L) |
| PGS001458 | Cingulum — cingulate gyrus (R) |
| PGS001459 | Cingulum — hippocampus (L) |
| PGS001460 | Cingulum — hippocampus (R) |
| PGS001466 | Genu of corpus callosum |
| PGS001471 | Middle cerebellar peduncle |
| PGS001474 | Posterior limb of internal capsule (L) |
| PGS001478 | Retrolenticular part of internal capsule (L) |
| PGS001479 | Retrolenticular part of internal capsule (R) |
| PGS001480 | Sagittal stratum (L) |
| PGS001481 | Sagittal stratum (R) |
| PGS001484 | Superior corona radiata (L) |
| PGS001485 | Superior corona radiata (R) |
| PGS001662 | Acoustic radiation (R) † |
| PGS001669 | Forceps major |
| PGS001679 | Parahippocampal part of cingulum (R) † |

† **Negative controls** — these PGS do not contain rs55705857 or its LD proxies in their scoring file. Conditioning on rs55705857 should *not* change their association.

## Extending to Molecular Subtypes

The pipeline supports three phenotype definitions via `--phenotype`:

| Phenotype | Cases | Controls | Biological question |
|-----------|-------|----------|-------------------|
| `idhmt` (default) | IDH-mutant glioma (idh=1) | Population controls (case=0) | All IDH-mutant gliomas |
| `idhmt_intact` | IDH-mutant, 1p/19q-intact (idh=1, pq=0) | Population controls | Astrocytoma specifically |
| `idhmt_codel` | IDH-mutant, 1p/19q-codeleted (idh=1, pq=1) | Population controls | Oligodendroglioma specifically |

Run all three with:

```bash
for pheno in idhmt idhmt_intact idhmt_codel; do
    bash run_pipeline.sh \
        --scores-dir /data/scores \
        --vcf-hg19-dir /data/imputed_hg19 \
        --output-dir ./output \
        --phenotype "$pheno" \
        --skip-extraction   # dosage is phenotype-independent
done
```

Note: `--skip-extraction` for the second and third runs because dosage extraction is the same across phenotypes.

## Pipeline Structure

```
project/
├── run_pipeline.sh              # Master runner (this file)
├── README.md                    # This documentation
├── scripts/
│   ├── 01_extract_dosage.sh     # Dosage extraction
│   ├── 02_merge_data.py         # Data merging
│   ├── 03_analysis1_stratified.py
│   ├── 04_analysis2_conditional.py
│   ├── 05_analysis3_interaction.py
│   ├── 06a_prep_ld_pruning.sh
│   ├── 06_analysis4_ld_pruned.py
│   ├── 07_figures.py
│   ├── 08_tables.py
│   └── 09_report.tex
└── utils/
    ├── __init__.py
    ├── params.py                # Shared constants (PGS IDs, variant info)
    └── analysis_utils.py        # Shared analysis functions
```

## Troubleshooting

### bcftools not found

Install via conda or your system package manager:
```bash
conda install -c bioconda bcftools
# or
sudo apt-get install bcftools
```
The pipeline falls back to `tabix + awk` if bcftools is unavailable, but bcftools is preferred for reliability.

### Variant not found at expected position

- Verify your VCFs use **hg19** coordinates and chromosome `8` (no `chr` prefix)
- Check that the VCF contains position 130645692 on chromosome 8:
  ```bash
  tabix your_vcf.gz 8:130645692-130645692
  ```
- If your VCFs use hg38, liftover to hg19 first or adjust the coordinate to `chr8:129633446`

### PGS column not found

- Check that PGS scoring was performed with all 18 ICVF models listed in `utils/params.py`
- The first column of the z-scores file must be named `sample` or `IID`
- The pipeline logs which PGS columns are missing vs available

### No samples after merge

- Check that sample IDs match across covariates, PGS, and dosage files
- Verify `exclude` column is not filtering everyone
- Check the `merge_summary.tsv` for per-step sample counts

### SLURM job fails

Check `output/{phenotype}/logs/slurm_*.out` for the error. Common issues:
- Insufficient memory: increase `--slurm-mem`
- Module not loaded: add `module load` commands to your environment
- Path issues: use absolute paths for all `--*-dir` arguments

### pdflatex errors

The LaTeX report step is optional. If compilation fails:
- Check `output/{phenotype}/logs/09_report_pass*.log`
- Ensure figure/table files exist in the expected locations
- Install missing LaTeX packages via `tlmgr install <package>`
