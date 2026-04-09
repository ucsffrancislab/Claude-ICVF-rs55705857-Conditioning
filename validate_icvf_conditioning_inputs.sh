#!/usr/bin/env bash
# =============================================================================
# validate_icvf_conditioning_inputs.sh
#
# Pre-flight checks for the ICVF conditioning analysis pipeline.
# Validates:
#   1. The 18 ICVF PGS exist in the scored data files
#   2. rs55705857 (chr8:128748020 on hg38) is present and well-imputed
#   3. PGS variant match rates (from scoring logs, if available)
#
# Usage:
#   bash validate_icvf_conditioning_inputs.sh \
#     --scores-dir /path/to/scores \
#     --vcf-dir    /path/to/hg38/vcfs \
#     [--vcf-hg19-dir /path/to/hg19/vcfs] \
#     [--raw-geno-dir /path/to/raw/genotypes] \
#     [--scoring-log-dir /path/to/scoring/logs] \
#     [--output validation_report.txt]
#
# The script is non-destructive (read-only) and prints a summary report.
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
SCORES_DIR=""
VCF_DIR=""
VCF_HG19_DIR=""
RAW_GENO_DIR=""
SCORING_LOG_DIR=""
OUTPUT="validation_report.txt"

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --scores-dir)      SCORES_DIR="$2";      shift 2 ;;
    --vcf-dir)         VCF_DIR="$2";         shift 2 ;;
    --vcf-hg19-dir)    VCF_HG19_DIR="$2";    shift 2 ;;
    --raw-geno-dir)    RAW_GENO_DIR="$2";    shift 2 ;;
    --scoring-log-dir) SCORING_LOG_DIR="$2"; shift 2 ;;
    --output)          OUTPUT="$2";          shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "$SCORES_DIR" || -z "$VCF_DIR" ]]; then
  echo "ERROR: --scores-dir and --vcf-dir are required."
  echo "Run with --help or see script header for usage."
  exit 1
fi

# ── Constants ─────────────────────────────────────────────────────────────────

# The 18 ICVF PGS IDs
ICVF_PGS=(
  PGS001454  # Mean ICVF in body of corpus callosum
  PGS001456  # Mean ICVF in cerebral peduncle (R)
  PGS001457  # Mean ICVF in cingulum cingulate gyrus (L)
  PGS001458  # Mean ICVF in cingulum cingulate gyrus (R)
  PGS001459  # Mean ICVF in cingulum hippocampus (L)
  PGS001460  # Mean ICVF in cingulum hippocampus (R)
  PGS001466  # Mean ICVF in genu of corpus callosum
  PGS001471  # Mean ICVF in middle cerebellar peduncle
  PGS001474  # Mean ICVF in posterior limb of internal capsule (L)
  PGS001478  # Mean ICVF in retrolenticular part of internal capsule (L)
  PGS001479  # Mean ICVF in retrolenticular part of internal capsule (R)
  PGS001480  # Mean ICVF in sagittal stratum (L)
  PGS001481  # Mean ICVF in sagittal stratum (R)
  PGS001484  # Mean ICVF in superior corona radiata (L)
  PGS001485  # Mean ICVF in superior corona radiata (R)
  PGS001662  # WA ICVF in tract acoustic radiation (R)
  PGS001669  # WA ICVF in tract forceps major
  PGS001679  # WA ICVF in tract parahippocampal part of cingulum (R)
)

# Expected variant counts per PGS (from PGS Catalog)
declare -A EXPECTED_VARIANTS
EXPECTED_VARIANTS[PGS001454]=2881
EXPECTED_VARIANTS[PGS001456]=1501
EXPECTED_VARIANTS[PGS001457]=1929
EXPECTED_VARIANTS[PGS001458]=3061
EXPECTED_VARIANTS[PGS001459]=1153
EXPECTED_VARIANTS[PGS001460]=2400
EXPECTED_VARIANTS[PGS001466]=2530
EXPECTED_VARIANTS[PGS001471]=1580
EXPECTED_VARIANTS[PGS001474]=2010
EXPECTED_VARIANTS[PGS001478]=2268
EXPECTED_VARIANTS[PGS001479]=2113
EXPECTED_VARIANTS[PGS001480]=1651
EXPECTED_VARIANTS[PGS001481]=2195
EXPECTED_VARIANTS[PGS001484]=2789
EXPECTED_VARIANTS[PGS001485]=2299
EXPECTED_VARIANTS[PGS001662]=979
EXPECTED_VARIANTS[PGS001669]=1466
EXPECTED_VARIANTS[PGS001679]=955

# rs55705857 coordinates
RS_HG38_CHR="chr8"      # or "8" — we check both
RS_HG38_POS="128748020"
RS_HG19_CHR="chr8"
RS_HG19_POS="130194364"  # hg19 coordinate for rs55705857
RS_ID="rs55705857"

# Dataset names (adjust if yours differ)
DATASETS=(cidr i370 onco tcga)

# ── Helper functions ──────────────────────────────────────────────────────────

PASS=0
WARN=0
FAIL=0

log()   { echo "$@" | tee -a "$OUTPUT"; }
pass()  { log "  ✅ PASS: $*"; ((PASS++)) || true; }
warn()  { log "  ⚠️  WARN: $*"; ((WARN++)) || true; }
fail()  { log "  ❌ FAIL: $*"; ((FAIL++)) || true; }
header(){ log ""; log "=============================================================================="; log "$*"; log "=============================================================================="; }

# ── Initialize report ─────────────────────────────────────────────────────────
> "$OUTPUT"
log "ICVF Conditioning Analysis — Input Validation Report"
log "Generated: $(date)"
log "Scores dir:       $SCORES_DIR"
log "VCF dir (hg38):   $VCF_DIR"
[[ -n "$VCF_HG19_DIR" ]]    && log "VCF dir (hg19):   $VCF_HG19_DIR"
[[ -n "$RAW_GENO_DIR" ]]    && log "Raw genotype dir:  $RAW_GENO_DIR"
[[ -n "$SCORING_LOG_DIR" ]] && log "Scoring log dir:   $SCORING_LOG_DIR"

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CHECK 1: Do the 18 ICVF PGS exist in the scored data?                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
header "CHECK 1: ICVF PGS presence in score files"

# Find score files — try common patterns
SCORE_FILES=()
for ds in "${DATASETS[@]}"; do
  for pattern in \
    "${SCORES_DIR}/${ds}.scores.z-scores.txt.gz" \
    "${SCORES_DIR}/${ds}.scores.z-scores.txt" \
    "${SCORES_DIR}/${ds}.sscore.zst" \
    "${SCORES_DIR}/${ds}_scores.txt.gz" \
    "${SCORES_DIR}/${ds}_scores.txt" \
    "${SCORES_DIR}/${ds}.scores.txt.gz" \
    "${SCORES_DIR}/${ds}.scores.txt" \
    "${SCORES_DIR}/${ds}/"*score* ; do
    if compgen -G "$pattern" > /dev/null 2>&1; then
      for f in $pattern; do
        SCORE_FILES+=("$f")
      done
      break
    fi
  done
done

if [[ ${#SCORE_FILES[@]} -eq 0 ]]; then
  fail "No score files found in $SCORES_DIR for datasets: ${DATASETS[*]}"
  log "  Tried patterns like {dataset}.scores.z-scores.txt.gz etc."
  log "  Please check --scores-dir path or adjust DATASETS array in this script."
else
  log "  Found ${#SCORE_FILES[@]} score file(s):"
  for f in "${SCORE_FILES[@]}"; do
    log "    $(basename "$f")  ($(du -h "$f" | cut -f1))"
  done

  # Extract header from first score file to check columns
  FIRST_SCORE="${SCORE_FILES[0]}"
  log ""
  log "  Checking columns in: $(basename "$FIRST_SCORE")"

  # Handle gzip
  if [[ "$FIRST_SCORE" == *.gz ]]; then
    HEADER=$(zcat "$FIRST_SCORE" | head -1)
  else
    HEADER=$(head -1 "$FIRST_SCORE")
  fi

  # Count total PGS columns (rough: anything with PGS in the name)
  TOTAL_PGS_COLS=$(echo "$HEADER" | tr '\t' '\n' | grep -ci 'PGS' || true)
  log "  Total PGS-related columns: $TOTAL_PGS_COLS"

  # Check each of the 18 ICVF PGS
  FOUND=0
  MISSING=()
  log ""
  log "  ICVF PGS column check:"
  for pgs in "${ICVF_PGS[@]}"; do
    # Case-insensitive grep in header; PGS IDs might appear as column names
    if echo "$HEADER" | tr '\t' '\n' | grep -qi "$pgs"; then
      ((FOUND++)) || true
    else
      MISSING+=("$pgs")
    fi
  done

  if [[ $FOUND -eq 18 ]]; then
    pass "All 18 ICVF PGS found in score file columns"
  elif [[ $FOUND -gt 0 ]]; then
    warn "Only $FOUND/18 ICVF PGS found. Missing: ${MISSING[*]}"
  else
    fail "None of the 18 ICVF PGS found in score file columns"
    log "  First 20 column names for reference:"
    echo "$HEADER" | tr '\t' '\n' | head -20 | while read -r col; do
      log "    $col"
    done
    log "    ..."
  fi

  # Show sample counts per file
  log ""
  log "  Sample counts per score file:"
  for f in "${SCORE_FILES[@]}"; do
    if [[ "$f" == *.gz ]]; then
      N=$(zcat "$f" | wc -l)
    else
      N=$(wc -l < "$f")
    fi
    log "    $(basename "$f"): $((N - 1)) samples"
  done
fi

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CHECK 2: rs55705857 in VCF files                                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
header "CHECK 2: rs55705857 in genotype data"

check_vcf_for_snp() {
  # $1 = label (e.g., "hg38")
  # $2 = directory
  # $3 = chr (e.g., "chr8" or "8")
  # $4 = position
  local label="$1" dir="$2" chr_prefix="$3" pos="$4"

  log ""
  log "  --- Searching $label VCFs in: $dir ---"

  local found_any=false

  for ds in "${DATASETS[@]}"; do
    # Try common VCF naming patterns for chr8
    local vcf=""
    for pattern in \
      "${dir}/${ds}"*chr8*.vcf.gz \
      "${dir}/${ds}"*chr8*.dose.vcf.gz \
      "${dir}/${ds}/"*chr8*.vcf.gz \
      "${dir}/chr8"*.vcf.gz \
      "${dir}/"*chr8*.dose.vcf.gz ; do
      if compgen -G "$pattern" > /dev/null 2>&1; then
        vcf=$(ls $pattern 2>/dev/null | head -1)
        break
      fi
    done

    if [[ -z "$vcf" ]]; then
      warn "$label / $ds: No chr8 VCF found (tried several patterns)"
      continue
    fi

    log "  $ds: checking $(basename "$vcf")"

    # Check if tabix index exists
    if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
      warn "$ds: No tabix/csi index for $(basename "$vcf") — falling back to grep (slower)"
      # Grep approach (slow but works)
      local result
      result=$(zgrep -m1 -E "^(chr)?8\s+${pos}\s" "$vcf" 2>/dev/null || true)
      if [[ -n "$result" ]]; then
        found_any=true
        # Extract key fields
        local ref alt qual info_field
        ref=$(echo "$result" | cut -f4)
        alt=$(echo "$result" | cut -f5)
        qual=$(echo "$result" | cut -f6)
        info_field=$(echo "$result" | cut -f8)
        # Try to extract R2/DR2/INFO
        local r2=""
        r2=$(echo "$info_field" | grep -oP '(R2|DR2|INFO)=[\d.]+' | head -1 || true)
        pass "$label / $ds: rs55705857 FOUND at pos $pos (REF=$ref ALT=$alt ${r2:-R2=?})"
      else
        fail "$label / $ds: rs55705857 NOT FOUND at pos $pos"
      fi
    else
      # Tabix approach (fast)
      # Try both "chr8" and "8" formats
      local result=""
      for chr_fmt in "chr8" "8"; do
        result=$(tabix "$vcf" "${chr_fmt}:${pos}-${pos}" 2>/dev/null || true)
        [[ -n "$result" ]] && break
      done

      if [[ -n "$result" ]]; then
        found_any=true
        local ref alt info_field r2 rsid
        rsid=$(echo "$result" | head -1 | cut -f3)
        ref=$(echo "$result" | head -1 | cut -f4)
        alt=$(echo "$result" | head -1 | cut -f5)
        info_field=$(echo "$result" | head -1 | cut -f8)
        r2=$(echo "$info_field" | grep -oP '(R2|DR2)=[\d.]+' | head -1 || true)
        pass "$label / $ds: rs55705857 FOUND (ID=$rsid REF=$ref ALT=$alt ${r2:-R2=?})"

        # Extract dosage stats (first 10 samples) to sanity-check allele freq
        local format_col gt_cols
        format_col=$(echo "$result" | head -1 | cut -f9)
        # Find DS field position in FORMAT
        local ds_idx=0 i=0
        IFS=':' read -ra fmt_fields <<< "$format_col"
        for field in "${fmt_fields[@]}"; do
          if [[ "$field" == "DS" ]]; then
            ds_idx=$i
            break
          fi
          ((i++)) || true
        done

        if [[ $ds_idx -gt 0 || "${fmt_fields[0]}" == "DS" ]]; then
          # Extract DS values from first 20 samples for allele frequency estimate
          local ds_values
          ds_values=$(echo "$result" | head -1 | cut -f10- | tr '\t' '\n' | head -50 | \
            cut -d: -f$((ds_idx+1)) | grep -v '^\.$' | head -20)
          if [[ -n "$ds_values" ]]; then
            local mean_ds
            mean_ds=$(echo "$ds_values" | awk '{s+=$1; n++} END {if(n>0) printf "%.4f", s/n; else print "NA"}')
            log "    Mean dosage (first 20 samples): $mean_ds (expected ~0.10-0.14 for RAF 5-7%)"
          fi
        fi
      else
        fail "$label / $ds: rs55705857 NOT FOUND at ${chr_fmt}:$pos"
      fi
    fi
  done

  if $found_any; then
    return 0
  else
    return 1
  fi
}

# Check hg38 VCFs (primary)
check_vcf_for_snp "hg38" "$VCF_DIR" "chr8" "$RS_HG38_POS"

# Check hg19 VCFs if provided
if [[ -n "$VCF_HG19_DIR" ]]; then
  check_vcf_for_snp "hg19" "$VCF_HG19_DIR" "chr8" "$RS_HG19_POS"
fi

# Check raw genotype data if provided (e.g., PLINK files)
if [[ -n "$RAW_GENO_DIR" ]]; then
  header "CHECK 2b: rs55705857 in raw (pre-imputation) genotype data"
  log "  Searching for $RS_ID in .bim/.pvar files..."

  found_raw=false
  for bim in "${RAW_GENO_DIR}/"*.bim "${RAW_GENO_DIR}/"*.pvar "${RAW_GENO_DIR}/"*/*.bim "${RAW_GENO_DIR}/"*/*.pvar; do
    [[ -f "$bim" ]] || continue
    if grep -q "$RS_ID\|128748020\|130194364" "$bim" 2>/dev/null; then
      pass "rs55705857 found in raw genotypes: $(basename "$bim")"
      grep "$RS_ID\|128748020\|130194364" "$bim" | head -3 | while read -r line; do
        log "    $line"
      done
      found_raw=true
      break
    fi
  done
  if ! $found_raw; then
    log "  Not found in raw genotype files — it was likely imputed (which is fine)."
  fi
fi

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CHECK 3: PGS variant match rates (from scoring logs)                    ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
header "CHECK 3: PGS scoring variant match rates"

if [[ -n "$SCORING_LOG_DIR" ]]; then
  log "  Searching scoring logs in: $SCORING_LOG_DIR"

  for pgs in "${ICVF_PGS[@]}"; do
    expected=${EXPECTED_VARIANTS[$pgs]}

    # Try to find a log file mentioning this PGS
    log_match=$(grep -rl "$pgs" "$SCORING_LOG_DIR" 2>/dev/null | head -1 || true)

    if [[ -n "$log_match" ]]; then
      # Try to extract variant counts — patterns vary by scoring software
      # plink2 --score outputs: "N variants processed", "N variants matched"
      matched=$(grep -A5 "$pgs" "$log_match" | grep -oiP '(\d+)\s*(variants?\s*)?match' | grep -oP '\d+' | head -1 || true)
      if [[ -n "$matched" ]]; then
        pct=$(awk "BEGIN {printf \"%.1f\", 100*$matched/$expected}")
        if (( $(echo "$matched/$expected > 0.80" | bc -l) )); then
          pass "$pgs: $matched/$expected variants matched (${pct}%)"
        elif (( $(echo "$matched/$expected > 0.50" | bc -l) )); then
          warn "$pgs: $matched/$expected variants matched (${pct}%) — moderate coverage"
        else
          fail "$pgs: $matched/$expected variants matched (${pct}%) — LOW coverage"
        fi
      else
        warn "$pgs: Found in logs but couldn't parse variant match count"
      fi
    else
      warn "$pgs: No scoring log found"
    fi
  done
else
  log "  No --scoring-log-dir provided. Skipping variant match rate check."
  log ""
  log "  ALTERNATIVE: If you used plink2 --score, check the .log files for lines like:"
  log "    --score: N variants processed, M variants matched"
  log ""
  log "  Or, if score files include a header with variant counts, we can check those."
  log ""
  log "  You can also verify match rates after the fact by checking the PGS scoring"
  log "  files (.txt.gz from PGS Catalog) against your imputed variant list."

  # Try a quick check: if score files have metadata or we can count non-NA values
  if [[ ${#SCORE_FILES[@]} -gt 0 ]]; then
    log ""
    log "  Quick sanity check — score value distributions (first score file):"
    FIRST_SCORE="${SCORE_FILES[0]}"

    for pgs in "${ICVF_PGS[@]}"; do
      if [[ "$FIRST_SCORE" == *.gz ]]; then
        HEADER_LINE=$(zcat "$FIRST_SCORE" | head -1)
      else
        HEADER_LINE=$(head -1 "$FIRST_SCORE")
      fi

      # Find column index for this PGS
      COL_IDX=$(echo "$HEADER_LINE" | tr '\t' '\n' | grep -n -i "$pgs" | head -1 | cut -d: -f1 || true)

      if [[ -n "$COL_IDX" ]]; then
        # Get basic stats (mean, sd, min, max, n_NA) from first 500 samples
        if [[ "$FIRST_SCORE" == *.gz ]]; then
          VALS=$(zcat "$FIRST_SCORE" | head -501 | tail -500 | cut -f"$COL_IDX")
        else
          VALS=$(head -501 "$FIRST_SCORE" | tail -500 | cut -f"$COL_IDX")
        fi
        STATS=$(echo "$VALS" | awk '
          BEGIN {n=0; s=0; ss=0; na=0; min=1e30; max=-1e30}
          {
            if ($1 == "NA" || $1 == "." || $1 == "") { na++; next }
            n++; v=$1+0; s+=v; ss+=v*v
            if (v<min) min=v; if (v>max) max=v
          }
          END {
            if (n>0) {
              mean=s/n; sd=sqrt(ss/n - mean*mean)
              printf "n=%d NA=%d mean=%.3f sd=%.3f range=[%.3f, %.3f]", n, na, mean, sd, min, max
            } else {
              printf "ALL NA or MISSING (%d rows)", na
            }
          }')
        log "    $pgs: $STATS"
      fi
    done
  fi
fi

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CHECK 4: Additional data inventory                                      ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
header "CHECK 4: Covariate / phenotype data inventory"

for ds in "${DATASETS[@]}"; do
  for pattern in \
    "${SCORES_DIR}/${ds}-covariates.csv" \
    "${SCORES_DIR}/${ds}_covariates.csv" \
    "${SCORES_DIR}/../${ds}-covariates.csv"; do
    if [[ -f "$pattern" ]]; then
      log "  $ds covariates: $pattern"
      # Check for IDH, age, sex, PCs
      COVAR_HEADER=$(head -1 "$pattern")
      for col in idh pq age sex PC1 PC8; do
        if echo "$COVAR_HEADER" | grep -qi "$col"; then
          log "    ✓ $col"
        else
          log "    ✗ $col — NOT FOUND"
        fi
      done
      # Count samples
      N=$(wc -l < "$pattern")
      log "    Samples: $((N - 1))"
      break
    fi
  done
done

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SUMMARY                                                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
header "SUMMARY"
log ""
log "  ✅ Passed:  $PASS"
log "  ⚠️  Warnings: $WARN"
log "  ❌ Failed:  $FAIL"
log ""

if [[ $FAIL -eq 0 && $WARN -eq 0 ]]; then
  log "  🎉 All checks passed! Ready to proceed with the ICVF conditioning pipeline."
elif [[ $FAIL -eq 0 ]]; then
  log "  ⚡ No critical failures, but review warnings above before proceeding."
else
  log "  🛑 Critical issues found. Please resolve failures before running the pipeline."
fi

log ""
log "  Report saved to: $OUTPUT"
log ""
log "  ─── NEXT STEPS ───"
log "  If all checks pass, share this report and I'll build the full pipeline."
log "  If PGS scores are missing, we'll need to score them using plink2 + PGS Catalog files."
log "  If rs55705857 is missing, we'll search for a high-LD proxy SNP."
