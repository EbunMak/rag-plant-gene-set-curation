#!/usr/bin/env bash
set -euo pipefail

################## GLOBAL VARIABLES #####################
THREADS=16
OUTDIR=$1    #"/scratch/42909494/gse118126_out"  #directory that fasta and meta downloads are in
METADIR="${OUTDIR}/meta"
FASTQDIR="${OUTDIR}/cleaned"
ALIGNDIR="${2}/align"   # new output directory on scratch for alignment "/scratch/jobid"
LOGDIR="log/"
STAR_INDEX=$3    #"/scratch/42930340/star_index/100/"  #this is where the star index folder is located (containing 'Genome' file)
STAR=STAR
SAMTOOLS=samtools
########################################################

echo "Alignment dir: ${ALIGNDIR}"

mkdir -p "${ALIGNDIR}" "${LOGDIR}"

SRR_LIST="${METADIR}/srr_list.txt"

[[ -f "${SRR_LIST}" ]] || { echo "Missing SRR list: ${SRR_LIST}"; exit 1; }
[[ -f "${STAR_INDEX}/Genome" ]] || { echo "Missing STAR index: ${STAR_INDEX}"; exit 1; }

echo "[STAR] Starting alignments..."

while read -r SRR; do

  R1=$(ls "${FASTQDIR}/${SRR}_1"*.fastq.gz)
  R2=$(ls "${FASTQDIR}/${SRR}_2"*.fastq.gz)

  [[ -f "${R1}" ]] || { echo "Missing ${R1}"; exit 1; }
  [[ -f "${R2}" ]] || { echo "Missing ${R2}"; exit 1; }

  PREFIX="${ALIGNDIR}/${SRR}."
  BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
  CSI="${BAM}.csi"

  if [[ -f "${BAM}" && -f "${CSI}" ]]; then
    echo "${SRR}: BAM and CSI exist, skipping."
    continue
  fi

  if [[ ! -f "${BAM}" ]]; then
    echo "Aligning ${SRR}..."

    ${STAR} \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${R1}" "${R2}" \
      --readFilesCommand zcat \
      --twopassMode Basic \
      --outFileNamePrefix "${PREFIX}" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM \
      --limitBAMsortRAM 350000000000  #350 GB
  else
    echo "${SRR}: BAM exists; making CSI index."
  fi

  ${SAMTOOLS} index -c "${BAM}"

done < "${SRR_LIST}"

echo "[STAR] Done."
