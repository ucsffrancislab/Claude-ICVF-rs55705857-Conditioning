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



