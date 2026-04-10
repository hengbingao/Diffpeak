#!/bin/bash --login
#SBATCH --job-name=diffpeak
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=hengbin.gao@uwa.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

### ====================================================
### diffpeak SLURM submission script
### Requires: pip install . (run once before submitting)
### Usage   : sbatch submit_diffpeak.sh
### ====================================================

### ====== Define paths ======

INPUT_DIR="/group/ll010/hgao/epimodifier/CUTnTag/singleepimod_comprehensive_CT/fragment_files"
OUTPUT_BASE="/group/ll010/hgao/epimodifier/CUTnTag/singleepimod_comprehensive_CT/diff_peaks"

### ====== Analysis parameters ======

LOG2FC=1.0        # Min |log2FC|
PVAL=0.05         # Max p-value    (replicates only)
QVAL=0.1          # Max q-value    (replicates only)
DIRECTION="both"  # both / up / down
GAP=0             # Peak merge gap (bp)
MIN_SIGNAL=0      # Min signal to include peak (filter background noise)
THREADS=8         # Match --cpus-per-task above

### ====== Activate environment ======

source activate
conda activate power

### ====== Run per histone mark per BEP ======

for histone in H2AK119ub1 H2BK20ac H3K14ac H3K18ac H3K27ac H3K27me1 H3K27me2 H3K27me3 \
               H3K36me2 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me2 H3K9me3 \
               H3R17me2a H3R8me2a H3S10ph H4K16ac H4K20ac H4K20me1 H4K20me3 \
               H4K5ac H4K8ac H4R3me2s Pol_II_Ser2P Total_RNA_Pol_II
do
    HISTONE_DIR="${INPUT_DIR}/${histone}"

    for BEP in BEP396_FOG1 BEP100_ZIM3 BEP137_p65HSF1 BEP486_SREBF2ddr
    do
        # ── Collect file lists ──────────────────────────────────────────
        TREAT_BG=$(ls    "${HISTONE_DIR}"/*.downsample.bgnorm 2>/dev/null | grep "${BEP}"      | tr '\n' ' ')
        CTRL_BG=$(ls     "${HISTONE_DIR}"/*.downsample.bgnorm 2>/dev/null | grep "BEP073_GFP"  | tr '\n' ' ')
        TREAT_PEAKS=$(ls "${HISTONE_DIR}"/*.downsample.bed    2>/dev/null | grep "${BEP}"      | tr '\n' ' ')
        CTRL_PEAKS=$(ls  "${HISTONE_DIR}"/*.downsample.bed    2>/dev/null | grep "BEP073_GFP"  | tr '\n' ' ')

        # ── Skip if any group is empty ──────────────────────────────────
        if [ -z "$TREAT_BG" ]; then
            echo "[SKIP] No treat bgnorm : ${histone} / ${BEP}"
            continue
        fi
        if [ -z "$CTRL_BG" ]; then
            echo "[SKIP] No ctrl bgnorm  : ${histone} / BEP073_GFP"
            continue
        fi
        if [ -z "$TREAT_PEAKS" ]; then
            echo "[SKIP] No treat peaks  : ${histone} / ${BEP}"
            continue
        fi
        if [ -z "$CTRL_PEAKS" ]; then
            echo "[SKIP] No ctrl peaks   : ${histone} / BEP073_GFP"
            continue
        fi

        # ── Output directory ────────────────────────────────────────────
        OUTDIR="${OUTPUT_BASE}/${histone}/${BEP}"
        mkdir -p "${OUTDIR}"

        echo "=============================="
        echo "[RUN] ${histone} / ${BEP}"

        # ── Run diffpeak (installed command) ────────────────────────────
        diffpeak \
            -t ${TREAT_BG} \
            -c ${CTRL_BG} \
            --treat-peaks ${TREAT_PEAKS} \
            --ctrl-peaks  ${CTRL_PEAKS} \
            --log2fc      ${LOG2FC} \
            --pval        ${PVAL} \
            --qval        ${QVAL} \
            --direction   ${DIRECTION} \
            --gap         ${GAP} \
            --min-signal  ${MIN_SIGNAL} \
            --threads     ${THREADS} \
            -o "${OUTDIR}" \
            -n "${histone}_${BEP}"

        if [ $? -eq 0 ]; then
            echo "[DONE] ${histone} / ${BEP}"
        else
            echo "[ERROR] diffpeak failed: ${histone} / ${BEP}" >&2
        fi

    done
done

echo "=============================="
echo "[diffpeak] All runs finished."
