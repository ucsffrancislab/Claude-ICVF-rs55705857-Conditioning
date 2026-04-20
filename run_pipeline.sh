#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master runner for the ICVF conditioning analysis pipeline
#
# Tests whether ICVF (brain white matter) polygenic scores become significant
# for IDH-mutant (or IDH-wildtype) glioma risk after conditioning on rs55705857.
#
# Usage:
#   bash run_pipeline.sh \
#       --scores-dir  /path/to/scores \
#       --vcf-hg19-dir /path/to/imputed_vcfs \
#       --output-dir  ./output \
#       --phenotype   idhmt
#
# See  bash run_pipeline.sh --help  for full options.
# =============================================================================
set -euo pipefail

# ── Locate pipeline root (where this script lives) ──────────────────────────
# Locate the directory this script lives in.
# Under SLURM, sbatch copies the script to a spool directory so dirname "$0"
# points to the wrong place. Use scontrol to recover the original path.
# Outside SLURM (interactive/local), dirname "$0" works fine.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Make the utils/ package importable from any working directory
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

# Load required modules (cluster-specific; safe no-op if module system absent)
if command -v module &>/dev/null; then
    module load plink2   2>/dev/null || true
    module load bcftools 2>/dev/null || true
    module load htslib   2>/dev/null || true   # provides tabix
fi

# ── Defaults ─────────────────────────────────────────────────────────────────
SCORES_DIR=""
VCF_HG19_DIR=""
LD_REF=""
PGS_CATALOG_DIR=""
OUTPUT_DIR="./output"
PHENOTYPE="idhmt"
DATASETS="cidr,i370,onco,tcga"
SKIP_EXTRACTION=0
SKIP_ANALYSES=0

# SLURM defaults
USE_SLURM=0
SLURM_PARTITION="normal"
SLURM_TIME="12:00:00"
SLURM_MEM="32G"
SLURM_CPUS=4

# ── Helpers ──────────────────────────────────────────────────────────────────
timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log_info()  { echo "[$(timestamp)] INFO  $*"; }
log_warn()  { echo "[$(timestamp)] WARN  $*" >&2; }
log_error() { echo "[$(timestamp)] ERROR $*" >&2; }

step_header() {
    local step_num="$1" step_name="$2"
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════╗"
    printf "║  Step %s: %-57s ║\n" "$step_num" "$step_name"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo ""
}

run_step() {
    # Usage: run_step <step_name> <log_file> <command> [args...]
    local step_name="$1"; shift
    local log_file="$1";  shift
    local start_time end_time elapsed

    start_time=$(date +%s)
    log_info "Running: $*"
    log_info "Log:     ${log_file}"

    if "$@" > >(tee -a "$log_file") 2> >(tee -a "$log_file" >&2); then
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        log_info "✓ ${step_name} completed in ${elapsed}s"
        return 0
    else
        local rc=$?
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        log_error "✗ ${step_name} FAILED (exit code ${rc}) after ${elapsed}s"
        log_error "  Check log: ${log_file}"
        return $rc
    fi
}

# ── Argument parsing ─────────────────────────────────────────────────────────
usage() {
    cat <<'EOF'
Usage: run_pipeline.sh [OPTIONS]

ICVF Conditioning Analysis Pipeline
====================================
Tests whether ICVF polygenic scores associate with glioma risk (IDH-mutant
or IDH-wildtype subtypes) after conditioning on the major risk variant
rs55705857 (OR ≈ 6.5).

Required arguments:
  --scores-dir DIR        Directory containing {dataset}-covariates.csv and
                          {dataset}.scores.z-scores.txt.gz files
  --vcf-hg19-dir DIR      Directory with per-dataset subdirs (cidr/, i370/,
                          onco/, tcga/) each containing chr8.dose.vcf.gz

Optional arguments:
  --ld-ref PATH           Path to plink2-format LD reference panel (e.g.
                          1000G EUR). Required for Analysis 4 (LD-pruned).
                          If omitted, Analysis 4 is skipped.
  --pgs-catalog-dir DIR   Path to PGS Catalog scoring files (.txt.gz).
                          Required only if --ld-ref is provided.
  --output-dir DIR        Top-level output directory (default: ./output)
  --phenotype STR         Phenotype: idhmt (default), idhmt_intact, idhmt_codel, idhwt, idhwt_gbm
  --datasets LIST         Comma-separated dataset names (default: cidr,i370,onco,tcga)
  --skip-extraction       Skip Step 1 (dosage extraction) if already done
  --skip-analyses         Skip Steps 3–6 and go straight to figures/tables

SLURM options:
  --slurm                 Submit as a SLURM batch job instead of running locally
  --slurm-partition STR   SLURM partition (default: normal)
  --slurm-time STR        SLURM wall time (default: 12:00:00)
  --slurm-mem STR         SLURM memory (default: 32G)
  --slurm-cpus INT        SLURM CPUs per task (default: 4)

  -h, --help              Show this help message
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scores-dir)       SCORES_DIR="$2";       shift 2 ;;
        --vcf-hg19-dir)     VCF_HG19_DIR="$2";     shift 2 ;;
        --ld-ref)           LD_REF="$2";            shift 2 ;;
        --pgs-catalog-dir)  PGS_CATALOG_DIR="$2";   shift 2 ;;
        --output-dir)       OUTPUT_DIR="$2";        shift 2 ;;
        --phenotype)        PHENOTYPE="$2";         shift 2 ;;
        --datasets)         DATASETS="$2";          shift 2 ;;
        --skip-extraction)  SKIP_EXTRACTION=1;      shift ;;
        --skip-analyses)    SKIP_ANALYSES=1;        shift ;;
        --slurm)            USE_SLURM=1;            shift ;;
        --slurm-partition)  SLURM_PARTITION="$2";   shift 2 ;;
        --slurm-time)       SLURM_TIME="$2";        shift 2 ;;
        --slurm-mem)        SLURM_MEM="$2";         shift 2 ;;
        --slurm-cpus)       SLURM_CPUS="$2";        shift 2 ;;
        -h|--help)          usage 0 ;;
        *)                  log_error "Unknown option: $1"; usage 1 ;;
    esac
done

# ── Validate required arguments ──────────────────────────────────────────────
if [[ -z "$SCORES_DIR" ]]; then
    log_error "--scores-dir is required"; usage 1
fi
if [[ -z "$VCF_HG19_DIR" ]] && [[ $SKIP_EXTRACTION -eq 0 ]]; then
    log_error "--vcf-hg19-dir is required (or use --skip-extraction)"; usage 1
fi
if [[ -n "$LD_REF" ]] && [[ -z "$PGS_CATALOG_DIR" ]]; then
    log_error "--pgs-catalog-dir is required when --ld-ref is provided"; usage 1
fi

# Validate phenotype
case "$PHENOTYPE" in
    idhmt|idhmt_intact|idhmt_codel|idhwt|idhwt_gbm) ;;
    *) log_error "Invalid phenotype: $PHENOTYPE (must be idhmt, idhmt_intact, idhmt_codel, idhwt, or idhwt_gbm)"; exit 1 ;;
esac

# Convert comma-separated datasets to space-separated for bash iteration
DATASETS_SPACE="${DATASETS//,/ }"

# ── Directory structure ──────────────────────────────────────────────────────
PHENO_DIR="${OUTPUT_DIR}/${PHENOTYPE}"
DOSAGE_DIR="${PHENO_DIR}/dosage"
MERGED_DIR="${PHENO_DIR}/merged"
RESULTS_DIR="${PHENO_DIR}/results"
FIGURES_DIR="${PHENO_DIR}/figures"
TABLES_DIR="${PHENO_DIR}/tables"
LOG_DIR="${PHENO_DIR}/logs"

mkdir -p "$DOSAGE_DIR" "$MERGED_DIR" "$RESULTS_DIR" "$FIGURES_DIR" \
         "$TABLES_DIR" "$LOG_DIR"

# ── Write config file for Python scripts ─────────────────────────────────────
CONFIG_FILE="${PHENO_DIR}/pipeline_config.tsv"
{
    echo -e "key\tvalue"
    echo -e "phenotype\t${PHENOTYPE}"
    echo -e "datasets\t${DATASETS}"
    echo -e "output_dir\t${OUTPUT_DIR}"
    echo -e "scores_dir\t${SCORES_DIR}"
    echo -e "vcf_hg19_dir\t${VCF_HG19_DIR}"
    echo -e "dosage_dir\t${DOSAGE_DIR}"
    echo -e "merged_dir\t${MERGED_DIR}"
    echo -e "results_dir\t${RESULTS_DIR}"
    echo -e "figures_dir\t${FIGURES_DIR}"
    echo -e "tables_dir\t${TABLES_DIR}"
    [[ -n "$LD_REF" ]] && echo -e "ld_ref\t${LD_REF}"
    [[ -n "$PGS_CATALOG_DIR" ]] && echo -e "pgs_catalog_dir\t${PGS_CATALOG_DIR}"
} > "$CONFIG_FILE"

# ── SLURM wrapper ────────────────────────────────────────────────────────────
if [[ $USE_SLURM -eq 1 ]]; then
    SLURM_SCRIPT="${LOG_DIR}/slurm_pipeline.sh"
    SLURM_LOG="${LOG_DIR}/slurm_%j.out"

    # Rebuild the command line without --slurm flags
    RERUN_ARGS=(
        "--scores-dir" "$SCORES_DIR"
        "--output-dir" "$OUTPUT_DIR"
        "--phenotype"  "$PHENOTYPE"
        "--datasets"   "$DATASETS"
    )
    [[ -n "$VCF_HG19_DIR" ]]    && RERUN_ARGS+=("--vcf-hg19-dir" "$VCF_HG19_DIR")
    [[ -n "$LD_REF" ]]          && RERUN_ARGS+=("--ld-ref" "$LD_REF")
    [[ -n "$PGS_CATALOG_DIR" ]] && RERUN_ARGS+=("--pgs-catalog-dir" "$PGS_CATALOG_DIR")
    [[ $SKIP_EXTRACTION -eq 1 ]] && RERUN_ARGS+=("--skip-extraction")
    [[ $SKIP_ANALYSES -eq 1 ]]   && RERUN_ARGS+=("--skip-analyses")

    cat > "$SLURM_SCRIPT" <<SLURM_EOF
#!/usr/bin/env bash
#SBATCH --job-name=icvf_pipeline_${PHENOTYPE}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --output=${SLURM_LOG}
#SBATCH --error=${SLURM_LOG}

set -euo pipefail
echo "SLURM Job ID: \${SLURM_JOB_ID}"
echo "Node: \$(hostname)"
echo "Start: \$(date)"

bash "${BASH_SOURCE[0]}" ${RERUN_ARGS[@]}

echo "End: \$(date)"
SLURM_EOF

    log_info "Submitting SLURM job..."
    log_info "Script: ${SLURM_SCRIPT}"
    sbatch "$SLURM_SCRIPT"
    exit 0
fi

# ── Print configuration ─────────────────────────────────────────────────────
PIPELINE_START=$(date +%s)

echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│           ICVF Conditioning Analysis Pipeline                   │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""
log_info "Configuration:"
log_info "  Scores dir:      ${SCORES_DIR}"
log_info "  VCF dir (hg19):  ${VCF_HG19_DIR:-<skipped>}"
log_info "  LD reference:    ${LD_REF:-<not provided — Analysis 4 will be skipped>}"
log_info "  PGS catalog dir: ${PGS_CATALOG_DIR:-<not provided>}"
log_info "  Output dir:      ${OUTPUT_DIR}"
log_info "  Phenotype:       ${PHENOTYPE}"
log_info "  Datasets:        ${DATASETS}"
log_info "  Config file:     ${CONFIG_FILE}"
log_info "  Skip extraction: ${SKIP_EXTRACTION}"
log_info "  Skip analyses:   ${SKIP_ANALYSES}"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1: Extract rs55705857 dosage
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_EXTRACTION -eq 0 ]]; then
    step_header "1/9" "Extract rs55705857 dosage from imputed VCFs"

    run_step "Dosage extraction" "${LOG_DIR}/01_extract_dosage.log" \
        bash "${SCRIPT_DIR}/scripts/01_extract_dosage.sh" \
            --vcf-dir  "$VCF_HG19_DIR" \
            --outdir   "$DOSAGE_DIR" \
            --datasets "$DATASETS_SPACE"
else
    log_info "Skipping Step 1 (dosage extraction) — using existing files in ${DOSAGE_DIR}"
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2: Merge covariates, PGS, and dosage
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_ANALYSES -eq 0 ]]; then
    step_header "2/9" "Merge data into analysis-ready files"

    run_step "Data merge" "${LOG_DIR}/02_merge_data.log" \
        python3 "${SCRIPT_DIR}/scripts/02_merge_data.py" \
            --scores-dir  "$SCORES_DIR" \
            --dosage-dir  "$DOSAGE_DIR" \
            --outdir      "$MERGED_DIR" \
            --phenotype   "$PHENOTYPE" \
            --datasets    $DATASETS_SPACE
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3: Analysis 1 — Stratified by rs55705857 carrier status
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_ANALYSES -eq 0 ]]; then
    step_header "3/9" "Analysis 1: Stratified by carrier status"

    run_step "Stratified analysis" "${LOG_DIR}/03_analysis1_stratified.log" \
        python3 "${SCRIPT_DIR}/scripts/03_analysis1_stratified.py" \
            --data-dir  "$MERGED_DIR" \
            --outdir      "$RESULTS_DIR" \
            --phenotype   "$PHENOTYPE"
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4: Analysis 2 — Conditional on rs55705857
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_ANALYSES -eq 0 ]]; then
    step_header "4/9" "Analysis 2: Conditional on rs55705857"

    run_step "Conditional analysis" "${LOG_DIR}/04_analysis2_conditional.log" \
        python3 "${SCRIPT_DIR}/scripts/04_analysis2_conditional.py" \
            --data-dir  "$MERGED_DIR" \
            --outdir      "$RESULTS_DIR" \
            --phenotype   "$PHENOTYPE"
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 5: Analysis 3 — PGS × rs55705857 interaction
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_ANALYSES -eq 0 ]]; then
    step_header "5/9" "Analysis 3: PGS × rs55705857 interaction"

    run_step "Interaction analysis" "${LOG_DIR}/05_analysis3_interaction.log" \
        python3 "${SCRIPT_DIR}/scripts/05_analysis3_interaction.py" \
            --data-dir  "$MERGED_DIR" \
            --outdir      "$RESULTS_DIR" \
            --phenotype   "$PHENOTYPE"
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 6: Analysis 4 — LD-pruned PGS (optional)
# ══════════════════════════════════════════════════════════════════════════════
if [[ $SKIP_ANALYSES -eq 0 ]]; then
    if [[ -n "$LD_REF" ]]; then
        step_header "6/9" "Analysis 4: LD-pruned PGS associations"

        # Step 6a: Prepare LD-pruned scoring files
        log_info "  6a: Preparing LD-pruned scoring files..."
        run_step "LD pruning prep" "${LOG_DIR}/06a_prep_ld_pruning.log" \
            bash "${SCRIPT_DIR}/scripts/06a_prep_ld_pruning.sh" \
                --ld-ref          "$LD_REF" \
                --pgs-catalog-dir "$PGS_CATALOG_DIR" \
                --vcf-hg19-dir   "$VCF_HG19_DIR" \
                --outdir          "$RESULTS_DIR/ld_pruned"

        # Step 6b: Run LD-pruned PGS analysis
        log_info "  6b: Running LD-pruned associations..."
        run_step "LD-pruned analysis" "${LOG_DIR}/06_analysis4_ld_pruned.log" \
            python3 "${SCRIPT_DIR}/scripts/06_analysis4_ld_pruned.py" \
                --data-dir      "$MERGED_DIR" \
                --pruned-scores-dir   "$RESULTS_DIR/ld_pruned" \
                --outdir          "$RESULTS_DIR" \
                --phenotype       "$PHENOTYPE"
    else
        step_header "6/9" "Analysis 4: LD-pruned PGS (SKIPPED)"
        log_warn "No --ld-ref provided. Skipping Analysis 4 (LD-pruned PGS)."
        log_warn "To run this analysis, provide --ld-ref and --pgs-catalog-dir."
    fi
fi

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 7: Publication figures
# ══════════════════════════════════════════════════════════════════════════════
step_header "7/9" "Generate publication figures"

run_step "Figures" "${LOG_DIR}/07_figures.log" \
    python3 "${SCRIPT_DIR}/scripts/07_figures.py" \
        --results-dir "$RESULTS_DIR" \
        --outdir      "$FIGURES_DIR" \
        --phenotype   "$PHENOTYPE"

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 8: Results tables
# ══════════════════════════════════════════════════════════════════════════════
step_header "8/9" "Generate results tables"

run_step "Tables" "${LOG_DIR}/08_tables.log" \
    python3 "${SCRIPT_DIR}/scripts/08_tables.py" \
        --data-dir    "$MERGED_DIR" \
        --results-dir "$RESULTS_DIR" \
        --outdir      "$TABLES_DIR" \
        --phenotype   "$PHENOTYPE"

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 9: Compile LaTeX report (optional)
# ══════════════════════════════════════════════════════════════════════════════
step_header "9/9" "Compile LaTeX report"

REPORT_TEX="${SCRIPT_DIR}/scripts/09_report.tex"
if [[ -f "$REPORT_TEX" ]] && command -v pdflatex &>/dev/null; then
    log_info "Compiling LaTeX report..."

    # Copy tex file to output and compile there
    cp "$REPORT_TEX" "${TABLES_DIR}/"
    REPORT_DIR="${TABLES_DIR}"
    (
        cd "$REPORT_DIR"
        # Run twice for cross-references
        pdflatex -interaction=nonstopmode "$(basename "$REPORT_TEX")" \
            > "${LOG_DIR}/09_report_pass1.log" 2>&1 || true
        pdflatex -interaction=nonstopmode "$(basename "$REPORT_TEX")" \
            > "${LOG_DIR}/09_report_pass2.log" 2>&1 || true
    )

    REPORT_PDF="${REPORT_DIR}/09_report.pdf"
    if [[ -f "$REPORT_PDF" ]]; then
        log_info "✓ Report compiled: ${REPORT_PDF}"
    else
        log_warn "LaTeX compilation did not produce a PDF."
        log_warn "Check ${LOG_DIR}/09_report_pass*.log for errors."
    fi
elif [[ ! -f "$REPORT_TEX" ]]; then
    log_warn "Report template not found: ${REPORT_TEX}"
    log_warn "Skipping report compilation."
else
    log_warn "pdflatex not found. Skipping report compilation."
    log_warn "Install texlive or run 'pdflatex 09_report.tex' manually."
fi

# ══════════════════════════════════════════════════════════════════════════════
#  Summary
# ══════════════════════════════════════════════════════════════════════════════
PIPELINE_END=$(date +%s)
PIPELINE_ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
PIPELINE_MIN=$(( PIPELINE_ELAPSED / 60 ))
PIPELINE_SEC=$(( PIPELINE_ELAPSED % 60 ))

echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│                    Pipeline Complete                            │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""
log_info "Total elapsed time: ${PIPELINE_MIN}m ${PIPELINE_SEC}s"
log_info ""
log_info "Output directory: ${PHENO_DIR}/"
log_info ""
log_info "Output structure:"
log_info "  ${PHENO_DIR}/"
log_info "  ├── dosage/           Extracted rs55705857 dosages per dataset"
log_info "  ├── merged/           Analysis-ready merged files per dataset"
log_info "  ├── results/          Analysis results (CSVs)"
log_info "  │   └── ld_pruned/    LD-pruned scoring files & results"
log_info "  ├── figures/          Publication-ready figures (PDF/PNG)"
log_info "  ├── tables/           Formatted results tables (CSV/LaTeX)"
log_info "  ├── logs/             Per-step log files"
log_info "  └── pipeline_config.tsv"
echo ""

# List key output files
log_info "Key output files:"
for f in \
    "${DOSAGE_DIR}/extraction_summary.tsv" \
    "${MERGED_DIR}/merge_summary.tsv" \
    "${RESULTS_DIR}/analysis1_stratified_results.csv" \
    "${RESULTS_DIR}/analysis2_conditional_results.csv" \
    "${RESULTS_DIR}/analysis3_interaction_results.csv" \
    "${RESULTS_DIR}/analysis4_ld_pruned_results.csv" \
    "${FIGURES_DIR}"/*.pdf \
    "${TABLES_DIR}"/*.csv \
; do
    if [[ -f "$f" ]]; then
        log_info "  ✓ $(basename "$f")  →  $f"
    fi
done
echo ""
log_info "Done."
