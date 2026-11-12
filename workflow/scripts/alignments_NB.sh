#!/usr/bin/env bash
set -euo pipefail

# ======


REF="/home/alessio/workspace/RAMD/resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa.mmi"

FASTQ="/home/alessio/workspace/RAMD/alignments/RUNX1_NB/RUNX1_NB.targeted_sub.fastq.gz"

OUTDIR="/home/alessio/workspace/RAMD/alignments/RUNX1_NB/targeted/roi"
LOGDIR="/home/alessio/workspace/RAMD/logs/RUNX1_NB"

SAMPLE="RUNX1_NB"

OUTBAM="${OUTDIR}/${SAMPLE}.minimap2.targeted.srt.bam"
LOGFILE="${LOGDIR}/${SAMPLE}.minimap2_align_targeted.log"

THREADS=16

# ======================
# DIRECTORY
# ======================

mkdir -p "$OUTDIR" "$LOGDIR"

echo "[`date`] Starting minimap2 alignment for ${SAMPLE}"

# ======================
# ALIGNMENT
# ======================

minimap2 -ax splice \
    -t ${THREADS} \
    ${REF} \
    ${FASTQ} 2> "${LOGFILE}" \
| samtools sort \
    -@ ${THREADS} \
    -O bam \
    -o "${OUTBAM}" 2>> "${LOGFILE}"

echo "[`date`] Alignment completed for ${SAMPLE}"

conda deactivate