#!/usr/bin/env bash
set -euo pipefail


#FASTQ_DIR="/scratch/42909494/gse118126_out/fastq"
#OUT_DIR="${FASTQ_DIR}/cleaned"
FASTQ_DIR=$1
OUT_DIR=$2

FASTP=fastp

mkdir -p ${OUT_DIR}
echo "Output Directory: ${OUT_DIR}"


## Confirm each fastq1 has a matching fastq2 for paired end reads ##
fastq1=($(ls ${FASTQ_DIR}/*_1.fastq.gz | sort))
fastq2=($(ls ${FASTQ_DIR}/*_2.fastq.gz | sort))

if [ ${#fastq1[@]} -eq ${#fastq2[@]} ]; then
    echo "OK: ${#fastq1[@]} R1 files and ${#fastq2[@]} R2 files"
else
    echo "ERROR: R1 count (${#fastq1[@]}) != R2 count (${#fastq2[@]})"
    exit 1
fi

#Print pairs to verify matching
for i in "${!fastq1[@]}"; do
    echo "Pair $((i+1)): ${fastq1[$i]} <-> ${fastq2[$i]}"
done


## Begin fastp cleaning
echo "[FASTP] Starting fastq cleaning..."

## Assumes PHRED33 scores
## adaptors autodetected
## default quality filtering parameters
for i in "${!fastq1[@]}"; do

    sample=$(basename ${fastq1[$i]} _1.fastq.gz)
    out1="${OUT_DIR}/${sample}_1_cleaned.fastq.gz"
    out2="${OUT_DIR}/${sample}_2_cleaned.fastq.gz"
    json="${OUT_DIR}/${sample}_fastp.json" 
    html="${OUT_DIR}/${sample}_fastp.html"

	if [[ -f "${out1}" && -f "${out2}" && -f "${json}" && -f "${html}" ]]; then
        echo "Files for ${sample} exist, skipping."
        continue
    fi

	echo "Processing $((i+1))/${#fastq1[@]}: ${sample}"

    ${FASTP} \
    -i ${fastq1[$i]} \
    -I ${fastq2[$i]} \
    -o "${out1}" \
    -O "${out2}" \
    --length_required 50 \
    --report_title ${sample} \
    --json "${json}"\
    --html "${html}" \
    --verbose

done

echo "Finished cleaning fastqs. Output is in ${OUT_DIR}"
