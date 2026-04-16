#!/usr/bin/env bash
# =============================================================================
# 01_extract_dosage.sh
#
# Extract rs55705857 dosage from per-dataset hg19 imputed VCFs.
#
# rs55705857 is at hg19 chr8:130645692 (REF=A, ALT=G, risk allele=G).
# VCFs use NO chr prefix (chromosome "8") and the ID column format is
# "8:130645692" rather than an rsID.
#
# For each dataset the script:
#   1. Locates chr8.dose.vcf.gz under <vcf-dir>/<dataset>/
#   2. Extracts the DS (dosage) field for the target variant
#   3. Validates REF=A, ALT=G
#   4. Reports allele frequency (AF) and imputation R² from INFO
#   5. Writes <outdir>/<dataset>_rs55705857_dosage.tsv  (IID \t rs55705857_dosage)
#
# Usage:
#   bash 01_extract_dosage.sh --vcf-dir /path/to/vcfs --outdir /path/to/out
#
# Dependencies: bcftools (preferred) OR tabix + awk
# =============================================================================
set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
VCF_DIR=""
OUTDIR=""
DATASETS="cidr i370 onco tcga"
CHROM="8"
POS=130645692
EXPECTED_REF="A"
EXPECTED_ALT="G"
VARIANT_LABEL="rs55705857"
REGION="${CHROM}:${POS}-${POS}"

# ── Logging helpers ──────────────────────────────────────────────────────────
timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log_info()  { echo "[$(timestamp)] INFO  $*"; }
log_warn()  { echo "[$(timestamp)] WARN  $*" >&2; }
log_error() { echo "[$(timestamp)] ERROR $*" >&2; }

# ── Argument parsing ─────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --vcf-dir DIR      Directory containing per-dataset subdirectories
                     (each with chr8.dose.vcf.gz and .tbi index)
  --outdir  DIR      Output directory for dosage TSV files

Optional:
  --datasets LIST    Space-separated dataset names (default: "cidr i370 onco tcga")
  -h, --help         Show this help message

Example:
  $(basename "$0") --vcf-dir /data/imputed --outdir /data/dosage
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf-dir)   VCF_DIR="$2";  shift 2 ;;
        --outdir)    OUTDIR="$2";   shift 2 ;;
        --datasets)  DATASETS="$2"; shift 2 ;;
        -h|--help)   usage 0 ;;
        *)           log_error "Unknown option: $1"; usage 1 ;;
    esac
done

if [[ -z "$VCF_DIR" ]]; then
    log_error "--vcf-dir is required"
    usage 1
fi
if [[ -z "$OUTDIR" ]]; then
    log_error "--outdir is required"
    usage 1
fi

mkdir -p "$OUTDIR"

# ── Tool detection ───────────────────────────────────────────────────────────
USE_BCFTOOLS=0
if command -v bcftools &>/dev/null; then
    USE_BCFTOOLS=1
    log_info "Using bcftools ($(bcftools --version | head -1))"
elif command -v tabix &>/dev/null; then
    log_info "bcftools not found; falling back to tabix + awk"
else
    log_error "Neither bcftools nor tabix found in PATH. Aborting."
    exit 1
fi

# ── Summary tracking ─────────────────────────────────────────────────────────
SUMMARY_FILE="${OUTDIR}/extraction_summary.tsv"
echo -e "dataset\tn_samples\tREF\tALT\tAF\tR2\tstatus" > "$SUMMARY_FILE"

FAIL_COUNT=0

# ── Per-dataset extraction ───────────────────────────────────────────────────
for DATASET in $DATASETS; do
    log_info "========== Processing dataset: ${DATASET} =========="

    VCF="${VCF_DIR}/${DATASET}/chr8.dose.vcf.gz"
    OUTFILE="${OUTDIR}/${DATASET}_${VARIANT_LABEL}_dosage.tsv"

    # --- Check input file ---
    if [[ ! -f "$VCF" ]]; then
        log_error "VCF not found: $VCF"
        echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tMISSING_VCF" >> "$SUMMARY_FILE"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # Check for index
    if [[ ! -f "${VCF}.tbi" ]] && [[ ! -f "${VCF%.gz}.csi" ]]; then
        log_warn "No .tbi or .csi index found for $VCF"
        log_warn "Attempting to index with tabix..."
        if command -v tabix &>/dev/null; then
            tabix -p vcf "$VCF"
        else
            log_error "Cannot index VCF without tabix. Skipping ${DATASET}."
            echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tNO_INDEX" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi
    fi

    # --- Extract variant metadata + dosages ---
    if [[ $USE_BCFTOOLS -eq 1 ]]; then
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # bcftools approach
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Step 1: Extract variant-level metadata (REF, ALT, AF, R2)
        # Try common INFO field names for AF and R² across imputation servers
        VARIANT_LINE=$(bcftools query \
            -r "$REGION" \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/R2\n' \
            "$VCF" 2>/dev/null || true)

        # If AF/R2 fields are missing, try alternative names (MAF, Rsq)
        if [[ -z "$VARIANT_LINE" ]]; then
            VARIANT_LINE=$(bcftools query \
                -r "$REGION" \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/Rsq\n' \
                "$VCF" 2>/dev/null || true)
        fi

        # Fallback: just get the basic fields
        if [[ -z "$VARIANT_LINE" ]]; then
            VARIANT_LINE=$(bcftools query \
                -r "$REGION" \
                -f '%CHROM\t%POS\t%REF\t%ALT\t.\t.\n' \
                "$VCF" 2>/dev/null || true)
        fi

        if [[ -z "$VARIANT_LINE" ]]; then
            log_error "Variant ${CHROM}:${POS} not found in ${VCF}"
            echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tVARIANT_NOT_FOUND" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi

        # Handle multiple records at the same position (pick the A>G one)
        MATCHED_LINE=""
        while IFS= read -r line; do
            LINE_REF=$(echo "$line" | cut -f3)
            LINE_ALT=$(echo "$line" | cut -f4)
            if [[ "$LINE_REF" == "$EXPECTED_REF" && "$LINE_ALT" == "$EXPECTED_ALT" ]]; then
                MATCHED_LINE="$line"
                break
            fi
        done <<< "$VARIANT_LINE"

        if [[ -z "$MATCHED_LINE" ]]; then
            log_error "Variant at ${CHROM}:${POS} found but alleles do not match expected REF=${EXPECTED_REF}, ALT=${EXPECTED_ALT}"
            log_error "Found: $(echo "$VARIANT_LINE" | head -1 | cut -f3-4)"
            echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tALLELE_MISMATCH" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi

        OBS_REF=$(echo "$MATCHED_LINE" | cut -f3)
        OBS_ALT=$(echo "$MATCHED_LINE" | cut -f4)
        OBS_AF=$(echo  "$MATCHED_LINE" | cut -f5)
        OBS_R2=$(echo  "$MATCHED_LINE" | cut -f6)

        log_info "Variant found: REF=${OBS_REF} ALT=${OBS_ALT} AF=${OBS_AF} R2=${OBS_R2}"

        # Step 2: Extract per-sample dosages
        # Header
        echo -e "IID\t${VARIANT_LABEL}_dosage" > "$OUTFILE"

        # Extract sample names and DS values
        bcftools query \
            -r "$REGION" \
            -i "REF=\"${EXPECTED_REF}\" && ALT=\"${EXPECTED_ALT}\"" \
            -f '[%SAMPLE\t%DS\n]' \
            "$VCF" >> "$OUTFILE"

        N_SAMPLES=$(( $(wc -l < "$OUTFILE") - 1 ))
        log_info "Extracted dosages for ${N_SAMPLES} samples → ${OUTFILE}"

        echo -e "${DATASET}\t${N_SAMPLES}\t${OBS_REF}\t${OBS_ALT}\t${OBS_AF}\t${OBS_R2}\tOK" >> "$SUMMARY_FILE"

    else
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # tabix + awk fallback approach
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Step 1: Pull the variant line(s)
        RAW_LINES=$(tabix "$VCF" "$REGION" 2>/dev/null || true)

        if [[ -z "$RAW_LINES" ]]; then
            log_error "Variant ${CHROM}:${POS} not found in ${VCF}"
            echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tVARIANT_NOT_FOUND" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi

        # Pick the line with REF=A, ALT=G
        MATCHED_LINE=$(echo "$RAW_LINES" | awk -F'\t' \
            -v ref="$EXPECTED_REF" -v alt="$EXPECTED_ALT" \
            '$4==ref && $5==alt {print; exit}')

        if [[ -z "$MATCHED_LINE" ]]; then
            log_error "Allele mismatch at ${CHROM}:${POS} in ${DATASET}"
            echo -e "${DATASET}\t0\tNA\tNA\tNA\tNA\tALLELE_MISMATCH" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi

        OBS_REF=$(echo "$MATCHED_LINE" | cut -f4)
        OBS_ALT=$(echo "$MATCHED_LINE" | cut -f5)

        # Parse AF and R2 from INFO field
        INFO_FIELD=$(echo "$MATCHED_LINE" | cut -f8)
        OBS_AF=$(echo "$INFO_FIELD" | grep -oP '(?:AF|MAF)=\K[^;]+' || echo ".")
        OBS_R2=$(echo "$INFO_FIELD" | grep -oP '(?:R2|Rsq)=\K[^;]+' || echo ".")

        log_info "Variant found: REF=${OBS_REF} ALT=${OBS_ALT} AF=${OBS_AF} R2=${OBS_R2}"

        # Step 2: Extract sample names from header
        HEADER_LINE=$(zgrep -m1 "^#CHROM" "$VCF")

        # Step 3: Find the DS index in FORMAT and extract dosages
        # Write header
        echo -e "IID\t${VARIANT_LABEL}_dosage" > "$OUTFILE"

        # Use awk to:
        #  - Determine DS field position from FORMAT column (col 9)
        #  - Extract that field from every sample column (col 10+)
        #  - Pair with sample names from the header
        paste <(echo "$HEADER_LINE") <(echo "$MATCHED_LINE") | awk -F'\t' '
        NR==1 {
            # First line is the header — store sample names
            for (i=10; i<=NF/2+0; i++) {
                # Header fields are in columns 1..NF/2
                # But we pasted header and data as one line, so header occupies
                # the first NF_hdr columns.
                # Actually let us do this differently below.
            }
        }
        ' 2>/dev/null  # placeholder — real extraction below

        # Simpler awk approach: process header and data separately
        # Get sample names (columns 10 onward from header)
        read -ra SAMPLES <<< "$(echo "$HEADER_LINE" | cut -f10-)"

        # Find DS index in FORMAT field
        FORMAT=$(echo "$MATCHED_LINE" | cut -f9)
        DS_INDEX=$(echo "$FORMAT" | tr ':' '\n' | grep -n "^DS$" | cut -d: -f1)

        if [[ -z "$DS_INDEX" ]]; then
            log_error "DS field not found in FORMAT column for ${DATASET}"
            echo -e "${DATASET}\t0\t${OBS_REF}\t${OBS_ALT}\t${OBS_AF}\t${OBS_R2}\tNO_DS_FIELD" >> "$SUMMARY_FILE"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            continue
        fi

        log_info "DS is field #${DS_INDEX} in FORMAT (${FORMAT})"

        # Extract dosage for each sample
        IDX=0
        for GENO in $(echo "$MATCHED_LINE" | cut -f10- | tr '\t' '\n'); do
            DOSAGE=$(echo "$GENO" | cut -d: -f"$DS_INDEX")
            echo -e "${SAMPLES[$IDX]}\t${DOSAGE}" >> "$OUTFILE"
            IDX=$((IDX + 1))
        done

        N_SAMPLES=$(( $(wc -l < "$OUTFILE") - 1 ))
        log_info "Extracted dosages for ${N_SAMPLES} samples → ${OUTFILE}"

        echo -e "${DATASET}\t${N_SAMPLES}\t${OBS_REF}\t${OBS_ALT}\t${OBS_AF}\t${OBS_R2}\tOK" >> "$SUMMARY_FILE"
    fi
done

# ── Final summary ────────────────────────────────────────────────────────────
log_info "========== Extraction Summary =========="
column -t -s$'\t' "$SUMMARY_FILE"

if [[ $FAIL_COUNT -gt 0 ]]; then
    log_warn "${FAIL_COUNT} dataset(s) failed. Check messages above."
    exit 1
fi

log_info "All datasets processed successfully."
log_info "Output directory: ${OUTDIR}"
exit 0
