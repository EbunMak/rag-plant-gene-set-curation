#!/usr/bin/env bash
set -euo pipefail

SRP="PRJNA950485"

THREADS=16

OUTDIR="prjna950485_out"

METADIR="${OUTDIR}/meta"
FASTQDIR="${OUTDIR}/fastq"

PYTHON=python3
PREFETCH=prefetch
FASTERQ=fasterq-dump

mkdir -p "${METADIR}" "${FASTQDIR}"

echo "[Download] Installing/checking pysradb..."
${PYTHON} -m pip install --user pysradb pandas

echo "[Download] Fetching metadata for ${SRP}..."

pysradb metadata "${SRP}" --desc --expand \
  > "${METADIR}/${SRP}.metadata.tsv"

echo "[Download] Extracting SRR run_accession list..."

awk -F'\t' '
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i=="run_accession") c=i
    }

    if (c=="") {
      print "ERROR: run_accession column not found" > "/dev/stderr"
      exit 1
    }
  }

  NR>1 {
    print $c
  }
' "${METADIR}/${SRP}.metadata.tsv" \
  > "${METADIR}/srr_list.txt"

echo "[Download] First few SRRs:"
head "${METADIR}/srr_list.txt"

echo "[Download] SRR count:"
wc -l "${METADIR}/srr_list.txt"

echo "[Download] Downloading SRA files..."

${PREFETCH} \
  --option-file "${METADIR}/srr_list.txt" \
  --output-directory "${FASTQDIR}"

echo "[Download] Converting to FASTQ..."

while read -r SRR; do

  if [[ -f "${FASTQDIR}/${SRR}_1.fastq.gz" && \
        -f "${FASTQDIR}/${SRR}_2.fastq.gz" ]]; then
    echo "${SRR}: FASTQs already exist, skipping."
    continue
  fi

  echo "======================================"
  echo "Processing ${SRR}"
  echo "======================================"

  ${FASTERQ} "${FASTQDIR}/${SRR}" \
    --split-files \
    --threads "${THREADS}" \
    --outdir "${FASTQDIR}"

  gzip -f "${FASTQDIR}/${SRR}_1.fastq"
  gzip -f "${FASTQDIR}/${SRR}_2.fastq"

  echo "${SRR}: complete."

done < "${METADIR}/srr_list.txt"

echo "[Download] Finished."
echo "FASTQs are in:"
echo "${FASTQDIR}"
