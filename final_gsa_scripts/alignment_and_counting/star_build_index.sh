#!/usr/bin/env bash
set -euo pipefail

####################### GLOBAL VARIABLES #########################
THREADS=16
#FASTA="ref/wheat.fa"
#GTF="ref/wheat.gtf"
FASTA="ref/Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.dna.toplevel.fa"
GTF="ref/Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.62.gtf"
SJDB_OVERHANG=$1    #read length -1  #exp. 149 or 100
SCRATCH=$2         #scratch directory # that files will be outputed to
STAR_INDEX="${SCRATCH}/star_index/${SJDB_OVERHANG}/"    #/scratch/slurmjobid/startIndex/overhang
echo ${STAR_INDEX}

STAR=STAR
###################################################################

echo "Read length -1:  ${SJDB_OVERHANG}"


mkdir -p "${STAR_INDEX}"

echo "[STAR] Checking reference files..."

[[ -f "${FASTA}" ]] || {
  echo "Missing FASTA: ${FASTA}"
  exit 1
}

[[ -f "${GTF}" ]] || {
  echo "Missing GTF: ${GTF}"
  exit 1
}

if [[ -f "${STAR_INDEX}/Genome" ]]; then
  echo "[STAR] Index already exists at ${STAR_INDEX}"
  exit 0
fi

echo "[STAR] Building genome index..."
echo "FASTA: ${FASTA}"
echo "GTF: ${GTF}"

${STAR} \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${FASTA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  --limitGenomeGenerateRAM 200000000000  #190GB   #190GB of RAM (it OOMed at 130GB)

echo "[STAR] Index build complete."
echo "Index location:"
echo "${STAR_INDEX}"
