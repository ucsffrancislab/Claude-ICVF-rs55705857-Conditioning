#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# 06a_prep_ld_pruning.sh — Prepare LD-pruned ICVF PGS scores.
#
# 1. Use plink2 to find variants in LD (r²>0.1) with rs55705857 (±1Mb).
# 2. For each ICVF PGS scoring file, identify and exclude overlapping variants.
# 3. Re-score each dataset with pruned scoring files via plink2 --score.
# 4. Standardise pruned scores to z-scores.
#
# Inputs:
#   --ld-ref          Plink2-format LD reference (prefix, e.g. 1000G EUR)
#   --pgs-catalog-dir Directory with PGS Catalog scoring files (.txt.gz)
#   --vcf-hg19-dir    Directory with per-dataset hg19 VCFs for re-scoring
#   --outdir          Output directory
#
# Outputs:
#   pruning_manifest.tsv   — per-PGS list of removed variant IDs
#   pruning_summary.tsv    — per-PGS removal counts
#   {dataset}_pruned_scores.tsv — per-dataset pruned PGS (z-scored)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Defaults ----
LD_REF=""
PGS_CATALOG_DIR=""
VCF_HG19_DIR=""
OUTDIR=""
LD_R2_THRESHOLD="0.1"
LD_WINDOW_KB="1000"
RS55705857_POS="130645692"       # hg19 position on chr8
RS55705857_ID="rs55705857"
CHR="8"

# ---- PGS IDs ----
PGS_IDS=(
  PGS001454 PGS001456 PGS001457 PGS001458 PGS001459 PGS001460
  PGS001466 PGS001471 PGS001474 PGS001478 PGS001479 PGS001480
  PGS001481 PGS001484 PGS001485 PGS001662 PGS001669 PGS001679
)

# PGS that do NOT contain rs55705857 (negative controls)
PGS_NO_RS=("PGS001456" "PGS001662" "PGS001679")

DATASETS=("cidr" "i370" "onco" "tcga")

# ---- Argument parsing ----
usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --ld-ref DIR            Plink2-format LD reference panel prefix
                          (e.g., path/to/1000G_EUR_chr8)
  --pgs-catalog-dir DIR   Directory with PGS Catalog scoring files (.txt.gz)
  --vcf-hg19-dir DIR      Directory with per-dataset hg19 VCFs
  --outdir DIR            Output directory

Optional:
  --ld-r2 FLOAT           LD r² threshold (default: 0.1)
  --ld-window-kb INT      LD window in kb (default: 1000)
  -h, --help              Show this message

If --ld-ref is not provided, the script prints instructions for obtaining
1000 Genomes EUR reference data and exits.
EOF
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ld-ref)          LD_REF="$2";          shift 2 ;;
    --pgs-catalog-dir) PGS_CATALOG_DIR="$2"; shift 2 ;;
    --vcf-hg19-dir)    VCF_HG19_DIR="$2";    shift 2 ;;
    --outdir)          OUTDIR="$2";           shift 2 ;;
    --ld-r2)           LD_R2_THRESHOLD="$2";  shift 2 ;;
    --ld-window-kb)    LD_WINDOW_KB="$2";     shift 2 ;;
    -h|--help)         usage 0 ;;
    *)                 echo "ERROR: Unknown option $1" >&2; usage 1 ;;
  esac
done

# ---- Logging helper ----
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }

# ---- Validate inputs ----
if [[ -z "$PGS_CATALOG_DIR" || -z "$VCF_HG19_DIR" || -z "$OUTDIR" ]]; then
  echo "ERROR: --pgs-catalog-dir, --vcf-hg19-dir, and --outdir are required." >&2
  usage 1
fi

if [[ -z "$LD_REF" ]]; then
  cat >&2 <<'INSTRUCTIONS'
========================================================================
  LD REFERENCE NOT PROVIDED
========================================================================
To run LD pruning you need a plink2-format LD reference for chr8.

Option A — 1000 Genomes Phase 3 EUR (recommended):
  1. Download:
     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
       ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
       ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

  2. Get EUR sample list:
     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
       integrated_call_samples_v3.20130502.ALL.panel
     awk '$3=="EUR" {print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel \
       > eur_samples.txt

  3. Convert to plink2 format (chr8 region around rs55705857):
     plink2 --vcf ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
       --keep eur_samples.txt \
       --chr 8 --from-bp 129645692 --to-bp 131645692 \
       --make-bed --out 1000G_EUR_chr8_rs55705857_region

  4. Then re-run this script with:
     --ld-ref 1000G_EUR_chr8_rs55705857_region

Option B — Use study controls:
  If you have plink2 files of your study genotypes, point --ld-ref to them.
  The script will extract the chr8 region automatically.
========================================================================
INSTRUCTIONS
  exit 1
fi

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/tmp_ld_pruning"
mkdir -p "${TMPDIR}"

log "Starting LD pruning pipeline"
log "  LD reference: ${LD_REF}"
log "  PGS catalog dir: ${PGS_CATALOG_DIR}"
log "  VCF dir: ${VCF_HG19_DIR}"
log "  Output dir: ${OUTDIR}"
log "  r² threshold: ${LD_R2_THRESHOLD}"
log "  LD window: ±${LD_WINDOW_KB}kb"

# =============================================================================
# STEP 1: Compute LD with rs55705857
# =============================================================================
log "STEP 1: Computing LD between rs55705857 and variants within ±${LD_WINDOW_KB}kb"

REGION_START=$(( RS55705857_POS - LD_WINDOW_KB * 1000 ))
REGION_END=$(( RS55705857_POS + LD_WINDOW_KB * 1000 ))

# Extract region from LD reference
plink2 \
  --bfile "${LD_REF}" \
  --chr "${CHR}" \
  --from-bp "${REGION_START}" \
  --to-bp "${REGION_END}" \
  --make-bed \
  --out "${TMPDIR}/region_extract" \
  2>&1 | tee "${TMPDIR}/plink2_extract.log" >&2

# Compute LD with rs55705857
# plink2 --ld-snp computes r² between the target SNP and all others
plink2 \
  --bfile "${TMPDIR}/region_extract" \
  --ld-snp "${RS55705857_ID}" \
  --ld-window-kb "${LD_WINDOW_KB}" \
  --ld-window-r2 0 \
  --r2-unphased \
  --out "${TMPDIR}/ld_rs55705857" \
  2>&1 | tee "${TMPDIR}/plink2_ld.log" >&2

# If plink2 --ld-snp fails (variant ID not found), try positional approach
if [[ ! -f "${TMPDIR}/ld_rs55705857.vcor2" ]]; then
  log "WARNING: --ld-snp by rsID failed; trying positional query"
  # Create a single-variant file for the target
  echo "${CHR}:${RS55705857_POS}" > "${TMPDIR}/target_var.txt"
  plink2 \
    --bfile "${TMPDIR}/region_extract" \
    --ld-snp-list "${TMPDIR}/target_var.txt" \
    --ld-window-kb "${LD_WINDOW_KB}" \
    --ld-window-r2 0 \
    --r2-unphased \
    --out "${TMPDIR}/ld_rs55705857" \
    2>&1 | tee -a "${TMPDIR}/plink2_ld.log" >&2
fi

# =============================================================================
# STEP 2: Extract variants exceeding r² threshold
# =============================================================================
log "STEP 2: Extracting variants with r² > ${LD_R2_THRESHOLD}"

LD_FILE="${TMPDIR}/ld_rs55705857.vcor2"
if [[ ! -f "$LD_FILE" ]]; then
  log "ERROR: LD output file not found: ${LD_FILE}"
  log "Check ${TMPDIR}/plink2_ld.log for details"
  exit 1
fi

# Parse the LD output — plink2 .vcor2 is tab-separated:
# columns: ID_A  ID_B  R2  (or similar)
# We want all ID_B where R2 > threshold (excluding the target itself)
LD_VARIANTS="${TMPDIR}/ld_variants_r2_above_threshold.txt"

python3 - "${LD_FILE}" "${LD_R2_THRESHOLD}" "${RS55705857_ID}" "${LD_VARIANTS}" <<'PYEOF'
import sys, csv, os

ld_file, r2_thresh, target_id, out_file = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]

variants = set()
# Try to parse plink2 LD output (tab-delimited, possibly with header)
with open(ld_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        # Typical plink2 .vcor2 columns: ID_A POS_A ID_B POS_B UNPHASED_R2
        # Or simpler: ID_A ID_B R2
        # We scan for the r2 value (last numeric column)
        try:
            r2 = float(parts[-1])
        except ValueError:
            continue
        if r2 >= r2_thresh:
            # Collect both IDs that are not the target
            for p in parts:
                if p != target_id:
                    try:
                        float(p)  # skip numeric values
                    except ValueError:
                        variants.add(p)

# Also store chr:pos pairs for positional matching
positions = set()
with open(ld_file) as f:
    for line in f:
        parts = line.strip().split()
        try:
            r2 = float(parts[-1])
        except ValueError:
            continue
        if r2 >= r2_thresh:
            for p in parts:
                if ':' in p:
                    positions.add(p)
                try:
                    pos_int = int(p)
                    if 128000000 < pos_int < 133000000:
                        positions.add(f"8:{pos_int}")
                except ValueError:
                    pass

with open(out_file, 'w') as f:
    for v in sorted(variants):
        f.write(v + '\n')

pos_file = out_file.replace('.txt', '_positions.txt')
with open(pos_file, 'w') as f:
    for p in sorted(positions):
        f.write(p + '\n')

print(f"Found {len(variants)} variant IDs with r2 >= {r2_thresh}", file=sys.stderr)
print(f"Found {len(positions)} chr:pos pairs", file=sys.stderr)
PYEOF

N_LD=$(wc -l < "${LD_VARIANTS}" || echo 0)
log "  Found ${N_LD} variants in LD (r² ≥ ${LD_R2_THRESHOLD}) with rs55705857"

# =============================================================================
# STEP 3: Prune PGS scoring files
# =============================================================================
log "STEP 3: Pruning PGS scoring files"

PRUNED_DIR="${OUTDIR}/pruned_scoring_files"
mkdir -p "${PRUNED_DIR}"

MANIFEST="${OUTDIR}/pruning_manifest.tsv"
SUMMARY="${OUTDIR}/pruning_summary.tsv"

# Header for manifest
echo -e "pgs_id\tvariant_id\tchr\tpos\teffect_allele\tweight" > "${MANIFEST}"
# Header for summary
echo -e "pgs_id\tn_total\tn_pruned\tn_remaining\tpct_removed\trs55705857_in_scoring" > "${SUMMARY}"

# Python helper to prune a single PGS scoring file
prune_pgs_file() {
  local PGS_ID="$1"
  local SCORING_FILE="$2"
  local LD_VAR_FILE="$3"
  local LD_POS_FILE="${LD_VAR_FILE%.txt}_positions.txt"
  local OUT_PRUNED="$4"
  local OUT_MANIFEST="$5"
  local OUT_SUMMARY="$6"

  python3 - "${PGS_ID}" "${SCORING_FILE}" "${LD_VAR_FILE}" "${LD_POS_FILE}" \
            "${OUT_PRUNED}" "${OUT_MANIFEST}" "${OUT_SUMMARY}" <<'PRUNE_PY'
import sys, gzip, os

pgs_id = sys.argv[1]
scoring_file = sys.argv[2]
ld_var_file = sys.argv[3]
ld_pos_file = sys.argv[4]
out_pruned = sys.argv[5]
out_manifest = sys.argv[6]
out_summary = sys.argv[7]

# Load LD variant IDs and positions
ld_ids = set()
if os.path.isfile(ld_var_file):
    with open(ld_var_file) as f:
        ld_ids = {line.strip() for line in f if line.strip()}

ld_positions = set()
if os.path.isfile(ld_pos_file):
    with open(ld_pos_file) as f:
        ld_positions = {line.strip() for line in f if line.strip()}

# Parse PGS Catalog scoring file
# Format: header lines starting with #, then tab-separated data
# Key columns: rsID (or chr_name:chr_position), effect_allele, effect_weight
opener = gzip.open if scoring_file.endswith('.gz') else open
header_line = None
data_lines = []
with opener(scoring_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if header_line is None:
            header_line = line.strip().split('	')
            continue
        data_lines.append(line.strip().split('	'))

# Find relevant column indices
col_map = {c.lower(): i for i, c in enumerate(header_line)}

# Possible ID columns
rsid_idx = col_map.get('rsid', col_map.get('snpid', col_map.get('variant_id', None)))
chr_idx = col_map.get('chr_name', col_map.get('chromosome', col_map.get('chr', None)))
pos_idx = col_map.get('chr_position', col_map.get('position', col_map.get('bp', col_map.get('pos', None))))
ea_idx = col_map.get('effect_allele', col_map.get('allele1', col_map.get('a1', None)))
weight_idx = col_map.get('effect_weight', col_map.get('weight', col_map.get('beta', None)))

n_total = len(data_lines)
pruned_variants = []
kept_lines = []

for row in data_lines:
    remove = False
    rsid_val = row[rsid_idx] if rsid_idx is not None and rsid_idx < len(row) else ''
    chr_val = row[chr_idx] if chr_idx is not None and chr_idx < len(row) else ''
    pos_val = row[pos_idx] if pos_idx is not None and pos_idx < len(row) else ''
    ea_val = row[ea_idx] if ea_idx is not None and ea_idx < len(row) else ''
    wt_val = row[weight_idx] if weight_idx is not None and weight_idx < len(row) else ''

    # Check if this variant should be removed
    if rsid_val in ld_ids:
        remove = True
    chr_pos = f"{chr_val}:{pos_val}"
    if chr_pos in ld_positions:
        remove = True
    # Also check without chr prefix
    chr_clean = chr_val.replace('chr', '')
    chr_pos_clean = f"{chr_clean}:{pos_val}"
    if chr_pos_clean in ld_positions:
        remove = True

    if remove:
        pruned_variants.append((rsid_val, chr_val, pos_val, ea_val, wt_val))
    else:
        kept_lines.append(row)

n_pruned = len(pruned_variants)
n_remaining = n_total - n_pruned
pct_removed = (n_pruned / n_total * 100) if n_total > 0 else 0.0

# Determine if rs55705857 was in the scoring file
rs_in = any(v[0] == 'rs55705857' for v in pruned_variants)
# Also check by position
if not rs_in:
    rs_in = any(v[2] == '130645692' for v in pruned_variants)

# Write pruned scoring file (same format, plink2 --score compatible)
with open(out_pruned, 'w') as f:
    f.write('	'.join(header_line) + '
')
    for row in kept_lines:
        f.write('	'.join(row) + '
')

# Append to manifest
with open(out_manifest, 'a') as f:
    for rsid_val, chr_val, pos_val, ea_val, wt_val in pruned_variants:
        f.write(f"{pgs_id}	{rsid_val}	{chr_val}	{pos_val}	{ea_val}	{wt_val}
")

# Append to summary
with open(out_summary, 'a') as f:
    rs_flag = 'yes' if rs_in else 'no'
    f.write(f"{pgs_id}	{n_total}	{n_pruned}	{n_remaining}	{pct_removed:.2f}	{rs_flag}
")

print(f"  {pgs_id}: {n_total} total, {n_pruned} pruned, {n_remaining} remaining ({pct_removed:.1f}% removed)", file=sys.stderr)
PRUNE_PY
}

for PGS_ID in "${PGS_IDS[@]}"; do
  # Find the scoring file (try common naming patterns)
  SCORING_FILE=""
  for pattern in "${PGS_CATALOG_DIR}/${PGS_ID}"*.txt.gz \
                  "${PGS_CATALOG_DIR}/${PGS_ID}"*.txt; do
    if [[ -f "$pattern" ]]; then
      SCORING_FILE="$pattern"
      break
    fi
  done

  if [[ -z "$SCORING_FILE" ]]; then
    log "WARNING: Scoring file not found for ${PGS_ID} in ${PGS_CATALOG_DIR} — skipping"
    echo -e "${PGS_ID}\t0\t0\t0\t0.00\tNA" >> "${SUMMARY}"
    continue
  fi

  PRUNED_FILE="${PRUNED_DIR}/${PGS_ID}_pruned.txt"

  prune_pgs_file "${PGS_ID}" "${SCORING_FILE}" "${LD_VARIANTS}" \
                 "${PRUNED_FILE}" "${MANIFEST}" "${SUMMARY}"
done

log "Pruning summary written to ${SUMMARY}"
log "Pruning manifest written to ${MANIFEST}"

# =============================================================================
# STEP 4: Re-score datasets with pruned scoring files
# =============================================================================
log "STEP 4: Re-scoring datasets with pruned PGS scoring files"

SCORES_DIR="${OUTDIR}/raw_pruned_scores"
mkdir -p "${SCORES_DIR}"

for DATASET in "${DATASETS[@]}"; do
  # Find VCF for this dataset
  VCF_FILE=""
  for pattern in "${VCF_HG19_DIR}/${DATASET}"*.vcf.gz \
                  "${VCF_HG19_DIR}/${DATASET}"*.vcf \
                  "${VCF_HG19_DIR}/${DATASET}"*.bcf; do
    if [[ -f "$pattern" ]]; then
      VCF_FILE="$pattern"
      break
    fi
  done

  if [[ -z "$VCF_FILE" ]]; then
    log "WARNING: VCF not found for ${DATASET} in ${VCF_HG19_DIR} — skipping"
    continue
  fi

  log "  Scoring dataset: ${DATASET} (${VCF_FILE})"

  for PGS_ID in "${PGS_IDS[@]}"; do
    PRUNED_FILE="${PRUNED_DIR}/${PGS_ID}_pruned.txt"
    if [[ ! -f "$PRUNED_FILE" ]]; then
      log "    Skipping ${PGS_ID} — no pruned scoring file"
      continue
    fi

    SCORE_PREFIX="${SCORES_DIR}/${DATASET}_${PGS_ID}"

    # Determine column indices for plink2 --score
    # PGS Catalog format: rsID chr_name chr_position effect_allele ...  effect_weight
    # We need: variant-ID-col  allele-col  score-col
    # Read header to figure out column numbers
    HEADER=$(head -1 "${PRUNED_FILE}")

    python3 -c "
import sys
header = '${HEADER}'.split('\t')
h = {c.lower(): i+1 for i, c in enumerate(header)}
# variant ID column
vid = h.get('rsid', h.get('snpid', h.get('variant_id', 1)))
# allele column
ea = h.get('effect_allele', h.get('allele1', h.get('a1', 4)))
# weight column
wt = h.get('effect_weight', h.get('weight', h.get('beta', len(header))))
print(f'{vid} {ea} {wt}')
" > "${TMPDIR}/col_indices_${PGS_ID}.txt"

    read -r VID_COL EA_COL WT_COL < "${TMPDIR}/col_indices_${PGS_ID}.txt"

    plink2 \
      --vcf "${VCF_FILE}" dosage=DS \
      --score "${PRUNED_FILE}" "${VID_COL}" "${EA_COL}" "${WT_COL}" \
        header cols=+scoresums \
      --out "${SCORE_PREFIX}" \
      2>&1 | tee "${SCORES_DIR}/${DATASET}_${PGS_ID}_score.log" >&2 || {
        log "    WARNING: plink2 --score failed for ${DATASET}/${PGS_ID}"
        continue
      }
  done
done

# =============================================================================
# STEP 5: Merge and z-score standardise pruned scores
# =============================================================================
log "STEP 5: Merging and standardising pruned scores"

for DATASET in "${DATASETS[@]}"; do
  python3 - "${DATASET}" "${SCORES_DIR}" "${OUTDIR}" <<'MERGE_PY'
import sys, os, glob
import numpy as np

dataset = sys.argv[1]
scores_dir = sys.argv[2]
outdir = sys.argv[3]

# Collect all .sscore files for this dataset
score_files = sorted(glob.glob(os.path.join(scores_dir, f"{dataset}_PGS*.sscore")))
if not score_files:
    print(f"  No score files found for {dataset} — skipping", file=sys.stderr)
    sys.exit(0)

# Read first file to get IIDs
merged = None
for sf in score_files:
    # Extract PGS ID from filename: {dataset}_{PGS_ID}.sscore
    basename = os.path.basename(sf)
    pgs_id = basename.replace(f"{dataset}_", "").replace(".sscore", "")

    # plink2 .sscore format: #IID  ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  SCORE1_AVG  SCORE1_SUM
    with open(sf) as f:
        header = f.readline().strip().split('\t')

    # Read as simple text parsing
    data = {}
    with open(sf) as f:
        hdr = f.readline().strip().split('\t')
        iid_idx = 0  # #IID or IID
        # Find score column (prefer SCORE1_SUM, else SCORE1_AVG)
        score_idx = None
        for i, h in enumerate(hdr):
            if 'SCORE1_SUM' in h.upper():
                score_idx = i
                break
        if score_idx is None:
            for i, h in enumerate(hdr):
                if 'SCORE1_AVG' in h.upper() or 'SCORE' in h.upper():
                    score_idx = i
                    break
        if score_idx is None:
            score_idx = len(hdr) - 1  # last column

        for line in f:
            parts = line.strip().split('\t')
            iid = parts[iid_idx].lstrip('#')
            try:
                val = float(parts[score_idx])
            except (ValueError, IndexError):
                val = float('nan')
            data[iid] = val

    col_name = f"{pgs_id}_pruned"
    if merged is None:
        merged = {iid: {'IID': iid, col_name: val} for iid, val in data.items()}
    else:
        for iid, val in data.items():
            if iid in merged:
                merged[iid][col_name] = val
            # Skip IIDs not in the first file

# Convert to list-of-dicts
rows = list(merged.values())
if not rows:
    print(f"  No data for {dataset}", file=sys.stderr)
    sys.exit(0)

# Get all columns
all_cols = sorted({k for r in rows for k in r if k != 'IID'})

# Z-score standardise each PGS column
for col in all_cols:
    vals = [r.get(col, float('nan')) for r in rows]
    arr = np.array(vals, dtype=float)
    valid = arr[np.isfinite(arr)]
    if len(valid) > 1:
        mu, sigma = np.mean(valid), np.std(valid, ddof=0)
        if sigma > 0:
            for r in rows:
                v = r.get(col, float('nan'))
                r[col] = (v - mu) / sigma if np.isfinite(v) else float('nan')

# Write output
out_path = os.path.join(outdir, f"{dataset}_pruned_scores.tsv")
with open(out_path, 'w') as f:
    header = ['IID'] + all_cols
    f.write('\t'.join(header) + '\n')
    for r in rows:
        vals = [r.get(c, '') for c in header]
        f.write('\t'.join(str(v) for v in vals) + '\n')

print(f"  {dataset}: {len(rows)} samples, {len(all_cols)} pruned PGS columns -> {out_path}", file=sys.stderr)
MERGE_PY
done

# =============================================================================
# Cleanup and summary
# =============================================================================
log "Cleaning up temporary files"
rm -rf "${TMPDIR}"

log "=== LD Pruning Pipeline Complete ==="
log "Outputs:"
log "  Pruning summary: ${SUMMARY}"
log "  Pruning manifest: ${MANIFEST}"
log "  Pruned scoring files: ${PRUNED_DIR}/"
for DATASET in "${DATASETS[@]}"; do
  OUTFILE="${OUTDIR}/${DATASET}_pruned_scores.tsv"
  if [[ -f "$OUTFILE" ]]; then
    log "  Pruned scores: ${OUTFILE}"
  fi
done

cat "${SUMMARY}"
