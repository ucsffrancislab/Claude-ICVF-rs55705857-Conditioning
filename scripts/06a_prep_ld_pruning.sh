#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# 06a_prep_ld_pruning.sh — Prepare LD-pruned ICVF PGS scores.
#
# 1. Use plink2 to find variants in LD (r²>0.1) with rs55705857 (±1Mb).
# 2. For each ICVF PGS scoring file, identify and exclude overlapping variants.
# 3. Re-score each dataset with pruned scoring files via plink2 --score,
#    looping over per-chromosome VCFs and aggregating SCORE1_SUM across chroms.
# 4. Standardise pruned scores to z-scores.
#
# Per-chromosome aggregation: a polygenic score is a plain weighted sum
#   PGS_i = sum_v w_v * d_{i,v}
# so it is mathematically identical to the sum of per-chromosome partial sums
#   PGS_i = sum_c sum_{v in c} w_v * d_{i,v}
# provided we aggregate SCORE1_SUM (the raw weighted sum), NOT SCORE1_AVG
# (which is (sum)/(2*nonmissing_v) and does not aggregate linearly).
#
# Inputs:
#   --ld-ref          Plink2-format LD reference (prefix, e.g. 1000G EUR)
#   --pgs-catalog-dir Directory with PGS Catalog scoring files (.txt.gz)
#   --vcf-hg19-dir    Directory with per-dataset hg19 per-chrom VCFs
#   --outdir          Output directory
#
# Outputs:
#   pruning_manifest.tsv         — per-PGS list of removed variant IDs
#   pruning_summary.tsv          — per-PGS removal counts
#   scoring_variant_coverage.tsv — per-dataset/PGS variant-match diagnostics
#   {dataset}_pruned_scores.tsv  — per-dataset aggregated pruned PGS (z-scored)
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

# Chromosomes to score (all autosomes)
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

# ---- Argument parsing ----
usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --ld-ref DIR            Plink2-format LD reference panel prefix
                          (e.g., path/to/1000G_EUR_chr8)
  --pgs-catalog-dir DIR   Directory with PGS Catalog scoring files (.txt.gz)
  --vcf-hg19-dir DIR      Directory with per-dataset per-chrom hg19 VCFs
                          (expected: {dir}/{dataset}/chr{N}.dose.vcf.gz)
  --outdir DIR            Output directory

Optional:
  --ld-r2 FLOAT           LD r² threshold (default: 0.1)
  --ld-window-kb INT      LD window in kb (default: 1000)
  -h, --help              Show this message
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
if [[ -z "$PGS_CATALOG_DIR" || -z "$VCF_HG19_DIR" || -z "$OUTDIR" || -z "$LD_REF" ]]; then
  echo "ERROR: --ld-ref, --pgs-catalog-dir, --vcf-hg19-dir, and --outdir are required." >&2
  usage 1
fi

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/tmp_ld_pruning"
mkdir -p "${TMPDIR}"

log "================================================================"
log "Starting LD pruning + per-chromosome PGS scoring pipeline"
log "  LD reference:    ${LD_REF}"
log "  PGS catalog dir: ${PGS_CATALOG_DIR}"
log "  VCF dir:         ${VCF_HG19_DIR}"
log "  Output dir:      ${OUTDIR}"
log "  r² threshold:    ${LD_R2_THRESHOLD}"
log "  LD window:       ±${LD_WINDOW_KB}kb"
log "  Chromosomes:     ${CHROMS[*]}"
log "  Datasets:        ${DATASETS[*]}"
log "  PGS count:       ${#PGS_IDS[@]}"
log "================================================================"

# =============================================================================
# STEP 1: Compute LD with rs55705857
# =============================================================================
log "STEP 1: Computing LD between rs55705857 and variants within ±${LD_WINDOW_KB}kb"

REGION_START=$(( RS55705857_POS - LD_WINDOW_KB * 1000 ))
REGION_END=$(( RS55705857_POS + LD_WINDOW_KB * 1000 ))

plink2 \
  --bfile "${LD_REF}" \
  --chr "${CHR}" \
  --from-bp "${REGION_START}" \
  --to-bp "${REGION_END}" \
  --make-bed \
  --out "${TMPDIR}/region_extract" \
  2>&1 | tee "${TMPDIR}/plink2_extract.log" >&2

plink2 \
  --bfile "${TMPDIR}/region_extract" \
  --ld-snp "${RS55705857_ID}" \
  --ld-window-kb "${LD_WINDOW_KB}" \
  --ld-window-r2 0 \
  --r2-unphased \
  --out "${TMPDIR}/ld_rs55705857" \
  2>&1 | tee "${TMPDIR}/plink2_ld.log" >&2

if [[ -f "${TMPDIR}/ld_rs55705857.vcor" ]]; then
  LD_FILE="${TMPDIR}/ld_rs55705857.vcor"
elif [[ -f "${TMPDIR}/ld_rs55705857.vcor2" ]]; then
  LD_FILE="${TMPDIR}/ld_rs55705857.vcor2"
else
  log "ERROR: LD output file not found (.vcor or .vcor2). See ${TMPDIR}/plink2_ld.log"
  exit 1
fi
log "  LD output: ${LD_FILE}"

# =============================================================================
# STEP 2: Extract variants exceeding r² threshold
# =============================================================================
log "STEP 2: Extracting variants with r² > ${LD_R2_THRESHOLD}"

LD_VARIANTS="${TMPDIR}/ld_variants_r2_above_threshold.txt"

python3 - "${LD_FILE}" "${LD_R2_THRESHOLD}" "${RS55705857_ID}" "${LD_VARIANTS}" <<'PYEOF'
import sys, os

ld_file, r2_thresh, target_id, out_file = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]

variants = set()
positions = set()
with open(ld_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        try:
            r2 = float(parts[-1])
        except ValueError:
            continue
        if r2 >= r2_thresh:
            for p in parts:
                if p == target_id:
                    continue
                try:
                    float(p)
                except ValueError:
                    variants.add(p)
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

print(f"  LD parse: {len(variants)} variant IDs, {len(positions)} chr:pos pairs", file=sys.stderr)
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

echo -e "pgs_id\tvariant_id\tchr\tpos\teffect_allele\tweight" > "${MANIFEST}"
echo -e "pgs_id\tn_total\tn_pruned\tn_remaining\tpct_removed\trs55705857_in_scoring" > "${SUMMARY}"

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

pgs_id, scoring_file, ld_var_file, ld_pos_file, out_pruned, out_manifest, out_summary = sys.argv[1:8]

ld_ids = set()
if os.path.isfile(ld_var_file):
    with open(ld_var_file) as f:
        ld_ids = {line.strip() for line in f if line.strip()}

ld_positions = set()
if os.path.isfile(ld_pos_file):
    with open(ld_pos_file) as f:
        ld_positions = {line.strip() for line in f if line.strip()}

opener = gzip.open if scoring_file.endswith('.gz') else open
header_line = None
data_lines = []
with opener(scoring_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if header_line is None:
            header_line = line.strip().split('\t')
            continue
        data_lines.append(line.strip().split('\t'))

col_map = {c.lower(): i for i, c in enumerate(header_line)}
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
    chr_val  = row[chr_idx]  if chr_idx  is not None and chr_idx  < len(row) else ''
    pos_val  = row[pos_idx]  if pos_idx  is not None and pos_idx  < len(row) else ''
    ea_val   = row[ea_idx]   if ea_idx   is not None and ea_idx   < len(row) else ''
    wt_val   = row[weight_idx] if weight_idx is not None and weight_idx < len(row) else ''

    if rsid_val in ld_ids:
        remove = True
    chr_pos = f"{chr_val}:{pos_val}"
    if chr_pos in ld_positions:
        remove = True
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

rs_in = any(v[0] == 'rs55705857' or v[2] == '130645692' for v in pruned_variants)

chr_col = header_line.index('chr_name') if 'chr_name' in header_line else 1
pos_col = header_line.index('chr_position') if 'chr_position' in header_line else 2

with open(out_pruned, 'w') as f:
    f.write('vcf_id\t' + '\t'.join(header_line) + '\n')
    for row in kept_lines:
        vcf_id = f"{row[chr_col]}:{row[pos_col]}"
        f.write(vcf_id + '\t' + '\t'.join(row) + '\n')

with open(out_manifest, 'a') as f:
    for rsid_val, chr_val, pos_val, ea_val, wt_val in pruned_variants:
        f.write(f"{pgs_id}\t{rsid_val}\t{chr_val}\t{pos_val}\t{ea_val}\t{wt_val}\n")

with open(out_summary, 'a') as f:
    rs_flag = 'yes' if rs_in else 'no'
    f.write(f"{pgs_id}\t{n_total}\t{n_pruned}\t{n_remaining}\t{pct_removed:.2f}\t{rs_flag}\n")

print(f"  {pgs_id}: {n_total} total, {n_pruned} pruned, {n_remaining} remaining ({pct_removed:.1f}% removed)", file=sys.stderr)
PRUNE_PY
}

for PGS_ID in "${PGS_IDS[@]}"; do
  SCORING_FILE=""
  for pattern in "${PGS_CATALOG_DIR}/${PGS_ID}"*.txt.gz \
                  "${PGS_CATALOG_DIR}/${PGS_ID}"*.txt; do
    if [[ -f "$pattern" ]]; then
      SCORING_FILE="$pattern"
      break
    fi
  done

  if [[ -z "$SCORING_FILE" ]]; then
    log "WARNING: Scoring file not found for ${PGS_ID} — skipping"
    echo -e "${PGS_ID}\t0\t0\t0\t0.00\tNA" >> "${SUMMARY}"
    continue
  fi

  PRUNED_FILE="${PRUNED_DIR}/${PGS_ID}_pruned.txt"
  prune_pgs_file "${PGS_ID}" "${SCORING_FILE}" "${LD_VARIANTS}" \
                 "${PRUNED_FILE}" "${MANIFEST}" "${SUMMARY}"
done

log "Pruning summary: ${SUMMARY}"
log "Pruning manifest: ${MANIFEST}"

# =============================================================================
# STEP 4: Re-score datasets per-chromosome, aggregate SCORE1_SUM across chroms
# =============================================================================
log "STEP 4: Per-chromosome scoring with pruned PGS files"
log "  --rm-dup exclude-all: drop all copies of duplicated variant IDs (e.g., multi-allelic"
log "                         sites collapsed to the same chr:pos). Principled & repeatable."
log "  Aggregation: sum of SCORE1_SUM (raw weighted sum) across 22 autosomes."
log "  SCORE1_AVG is NOT summed — it is a per-allele average whose denominator"
log "  differs per chrom and cannot be linearly combined."

SCORES_DIR="${OUTDIR}/raw_pruned_scores"
mkdir -p "${SCORES_DIR}"

# Diagnostic: per-dataset/PGS/chrom variant match counts
COVERAGE_FILE="${OUTDIR}/scoring_variant_coverage.tsv"
echo -e "dataset\tpgs_id\tchrom\tvariants_processed\tvariants_skipped\tn_dups_removed\tstatus" > "${COVERAGE_FILE}"

for DATASET in "${DATASETS[@]}"; do
  log "  ============================================"
  log "  Dataset: ${DATASET}"
  log "  ============================================"

  # Sanity-check that all expected per-chrom VCFs exist
  DS_VCF_DIR="${VCF_HG19_DIR}/${DATASET}"
  if [[ ! -d "$DS_VCF_DIR" ]]; then
    log "  WARNING: VCF dir ${DS_VCF_DIR} not found — skipping ${DATASET}"
    continue
  fi

  MISSING_CHROMS=()
  for CHR_N in "${CHROMS[@]}"; do
    if [[ ! -f "${DS_VCF_DIR}/chr${CHR_N}.dose.vcf.gz" ]]; then
      MISSING_CHROMS+=("${CHR_N}")
    fi
  done
  if [[ ${#MISSING_CHROMS[@]} -gt 0 ]]; then
    log "    NOTE: ${DATASET} missing chroms: ${MISSING_CHROMS[*]}"
  fi

  for PGS_ID in "${PGS_IDS[@]}"; do
    PRUNED_FILE="${PRUNED_DIR}/${PGS_ID}_pruned.txt"
    if [[ ! -f "$PRUNED_FILE" ]]; then
      log "    Skipping ${PGS_ID} — no pruned scoring file"
      continue
    fi

    # Determine column indices (1-based) for plink2 --score
    HEADER=$(head -1 "${PRUNED_FILE}")
    read -r VID_COL EA_COL WT_COL < <(python3 -c "
import sys
header = '''${HEADER}'''.split('\t')
h = {c.lower(): i+1 for i, c in enumerate(header)}
vid = h.get('vcf_id', 1)
ea  = h.get('effect_allele', 5)
wt  = h.get('effect_weight', 7)
print(f'{vid} {ea} {wt}')
")

    TOTAL_PROCESSED=0
    TOTAL_SKIPPED=0
    CHROMS_SCORED=0

    # --- Per-chromosome scoring ---
    for CHR_N in "${CHROMS[@]}"; do
      VCF_FILE="${DS_VCF_DIR}/chr${CHR_N}.dose.vcf.gz"
      if [[ ! -f "$VCF_FILE" ]]; then
        echo -e "${DATASET}\t${PGS_ID}\t${CHR_N}\t0\t0\tvcf_missing" >> "${COVERAGE_FILE}"
        continue
      fi

      SCORE_PREFIX="${SCORES_DIR}/${DATASET}_${PGS_ID}_chr${CHR_N}"
      SCORE_LOG="${SCORES_DIR}/${DATASET}_${PGS_ID}_chr${CHR_N}_score.log"

      plink2 \
        --vcf "${VCF_FILE}" dosage=DS \
        --rm-dup exclude-all list \
        --score "${PRUNED_FILE}" "${VID_COL}" "${EA_COL}" "${WT_COL}" \
          header cols=+scoresums \
        --out "${SCORE_PREFIX}" \
        > "${SCORE_LOG}" 2>&1 || {
          echo -e "${DATASET}\t${PGS_ID}\t${CHR_N}\t0\t0\t0\tplink2_failed" >> "${COVERAGE_FILE}"
          continue
        }

      # Parse variant counts from score log
      N_PROCESSED=$(grep -oP '\-\-score: \K[0-9]+(?= variants processed)' "${SCORE_LOG}" | tail -1)
      N_SKIPPED=$(grep -oP '\K[0-9]+(?= entries in)' "${SCORE_LOG}" | tail -1)
      # --rm-dup writes a "--rm-dup: N duplicate..." line and a .rmdup.list file
      N_DUPS=$(grep -oP '\-\-rm-dup: \K[0-9]+(?= duplicate)' "${SCORE_LOG}" | tail -1)
      N_PROCESSED="${N_PROCESSED:-0}"
      N_SKIPPED="${N_SKIPPED:-0}"
      N_DUPS="${N_DUPS:-0}"
      echo -e "${DATASET}\t${PGS_ID}\t${CHR_N}\t${N_PROCESSED}\t${N_SKIPPED}\t${N_DUPS}\tok" >> "${COVERAGE_FILE}"

      TOTAL_PROCESSED=$(( TOTAL_PROCESSED + N_PROCESSED ))
      CHROMS_SCORED=$(( CHROMS_SCORED + 1 ))
    done

    log "    ${PGS_ID}: ${TOTAL_PROCESSED} variants scored across ${CHROMS_SCORED} chromosomes"
  done
done

log "Per-chromosome scoring complete."
log "Variant coverage diagnostics: ${COVERAGE_FILE}"

# =============================================================================
# STEP 5: Aggregate SCORE1_SUM across chroms, then z-score standardise
# =============================================================================
log "STEP 5: Aggregating per-chrom SCORE1_SUM and standardising"

for DATASET in "${DATASETS[@]}"; do
  python3 - "${DATASET}" "${SCORES_DIR}" "${OUTDIR}" <<'MERGE_PY'
import sys, os, glob, re
import numpy as np

dataset, scores_dir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]

# Match per-chrom sscore files: {dataset}_{PGS_ID}_chr{N}.sscore
pattern = os.path.join(scores_dir, f"{dataset}_PGS*_chr*.sscore")
sscore_files = sorted(glob.glob(pattern))
if not sscore_files:
    print(f"  {dataset}: no sscore files found — skipping", file=sys.stderr)
    sys.exit(0)

# Group by PGS_ID
fname_re = re.compile(rf"{re.escape(dataset)}_(PGS\d+)_chr(\d+)\.sscore$")
by_pgs = {}
for sf in sscore_files:
    m = fname_re.search(os.path.basename(sf))
    if not m:
        continue
    pgs_id, chrom = m.group(1), int(m.group(2))
    by_pgs.setdefault(pgs_id, []).append((chrom, sf))

def read_sscore(path):
    """Return dict {iid: score_sum} using SCORE1_SUM column."""
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        # Locate IID column (first field, may be '#IID' or 'IID')
        iid_idx = 0
        # Locate SCORE1_SUM (required — we refuse to aggregate AVG)
        sum_idx = None
        for i, h in enumerate(hdr):
            if h.upper() == 'SCORE1_SUM':
                sum_idx = i
                break
        if sum_idx is None:
            # Fallback: look for any *_SUM column
            for i, h in enumerate(hdr):
                if h.upper().endswith('_SUM') and 'SCORE' in h.upper():
                    sum_idx = i
                    break
        if sum_idx is None:
            raise RuntimeError(f"No SCORE1_SUM column in {path} — cannot aggregate")

        data = {}
        for line in f:
            parts = line.rstrip('\n').split('\t')
            iid = parts[iid_idx].lstrip('#')
            try:
                data[iid] = float(parts[sum_idx])
            except (ValueError, IndexError):
                data[iid] = float('nan')
    return data

# Aggregate: per-PGS, sum SCORE1_SUM across chroms per IID
# This is mathematically identical to scoring against a single merged VCF,
# because PGS = sum_v w_v * d_{i,v} is additive over disjoint variant sets.
merged_rows = {}  # iid -> {col: value}
pgs_col_order = []

n_pgs_expected = len(by_pgs)
print(f"  {dataset}: aggregating {n_pgs_expected} PGS across per-chrom sscore files", file=sys.stderr)

for pgs_id in sorted(by_pgs.keys()):
    chrom_files = sorted(by_pgs[pgs_id])
    chrom_dicts = []
    for chrom, sf in chrom_files:
        chrom_dicts.append((chrom, read_sscore(sf)))

    # Sum across chroms — union of IIDs (should be identical across chroms)
    all_iids = set()
    for _, d in chrom_dicts:
        all_iids.update(d.keys())

    col_name = f"{pgs_id}_pruned"
    pgs_col_order.append(col_name)

    for iid in all_iids:
        total = 0.0
        contributing = 0
        for _, d in chrom_dicts:
            v = d.get(iid)
            if v is not None and np.isfinite(v):
                total += v
                contributing += 1
        merged_rows.setdefault(iid, {'IID': iid})[col_name] = (
            total if contributing > 0 else float('nan')
        )

    print(f"    {pgs_id}: aggregated {len(chrom_dicts)} chroms, "
          f"{len(all_iids)} samples", file=sys.stderr)

rows = list(merged_rows.values())
if not rows:
    print(f"  {dataset}: no data after aggregation", file=sys.stderr)
    sys.exit(0)

# Z-score standardise each PGS column (across samples)
for col in pgs_col_order:
    vals = np.array([r.get(col, float('nan')) for r in rows], dtype=float)
    valid = vals[np.isfinite(vals)]
    if len(valid) > 1:
        mu, sigma = float(np.mean(valid)), float(np.std(valid, ddof=0))
        if sigma > 0:
            for r in rows:
                v = r.get(col, float('nan'))
                r[col] = (v - mu) / sigma if np.isfinite(v) else float('nan')

out_path = os.path.join(outdir, f"{dataset}_pruned_scores.tsv")
with open(out_path, 'w') as f:
    header = ['IID'] + pgs_col_order
    f.write('\t'.join(header) + '\n')
    for r in rows:
        f.write('\t'.join(str(r.get(c, '')) for c in header) + '\n')

print(f"  {dataset}: wrote {len(rows)} samples × {len(pgs_col_order)} PGS -> {out_path}",
      file=sys.stderr)
MERGE_PY
done

# =============================================================================
# Cleanup and summary
# =============================================================================
log "STEP 6: Cleanup"
rm -rf "${TMPDIR}"

log "================================================================"
log "=== LD Pruning + Per-Chrom Scoring Pipeline Complete ==="
log "  Pruning summary:    ${SUMMARY}"
log "  Pruning manifest:   ${MANIFEST}"
log "  Pruned scoring dir: ${PRUNED_DIR}/"
log "  Per-chrom sscores:  ${SCORES_DIR}/"
log "  Coverage diag:      ${COVERAGE_FILE}"
for DATASET in "${DATASETS[@]}"; do
  OUTFILE="${OUTDIR}/${DATASET}_pruned_scores.tsv"
  if [[ -f "$OUTFILE" ]]; then
    log "  Aggregated scores:  ${OUTFILE}"
  fi
done
log "================================================================"

# Echo pruning summary to stdout for quick review
cat "${SUMMARY}"
