#!/usr/bin/env bash
set -euo pipefail

################# GLOBAL VARIABLES ####################################
THREADS=16

#OUTDIR="prjna950485_out"  #
#BAMDIR="${OUTDIR}/align"
OUTDIR=$1
BAMDIR=$2
COUNTSDIR="${OUTDIR}/counts"
GTFFILE="ref/Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.62.gtf"  #refseq v2.1
FEATURECOUNTS=featureCounts
######################################################################


# Create directories and ensure files exist
mkdir -p "${COUNTSDIR}"
echo "Created Counts Directory: ${COUNTSDIR}"

echo "BAM FILES:"
ls ${BAMDIR}/SRR*.Aligned.sortedByCoord.out.bam


[[ -f "${GTFFILE}" ]] || { echo "Missing GTFFILE: ${GTFFILE}"; exit 1; }
[[ -d "${BAMDIR}"  ]] || { echo "Error BAM file Directory: ${BAMDIR} not found"; exit 1; }


## RUN FEATURE COUNTS
echo "[FEATURECOUNTS] Starting feature counts..."

## THis feature count is for paired end reads and unstranded DNA, only counts reads that have both paired ends aligned
## It aslo counts multimapping reads (fractionally) which is beneficial for the hexaploid wheat genome. (will be rounded for DESEQ2 analysis)
${FEATURECOUNTS} \
    -T ${THREADS} \
    -a ${GTFFILE} \
    -o "${COUNTSDIR}/read_counts.txt" \
    -p \
    --countReadPairs \
    -B \
    -s 0 \
    -M \
    --fraction \
    --verbose \
    ${BAMDIR}/SRR*.Aligned.sortedByCoord.out.bam


echo "[FEATURECOUNTS] Done"

