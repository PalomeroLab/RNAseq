#!/bin/sh
## Template script for aligning reads with hisat2

SAMPLE_ID="$1"
INPUT_FOLDER="$2"
OUTPUT_FOLDER="$3"
THREADS="$4"
REFERENCE="$5"

# Derive the R1 and R2 FASTQ file names from sampleID
r1_fastq="${INPUT_FOLDER}/${SAMPLE_ID}_R1*.fastq.gz"
r2_fastq="${INPUT_FOLDER}/${SAMPLE_ID}_R2*.fastq.gz"

if [ ! -f "$r1_fastq" ] || [ ! -f "$r2_fastq" ]; then
	echo "Error: Missing R1 or R2 file for sample $SAMPLE_ID"
	return 1
fi

# Create output BAM file name
bam_file="${OUTPUT_FOLDER}/${SAMPLE_ID}_sorted_markdup.bam"

# Align reads and process BAM file
# the `dta` flag is only if we're aligning RNAseq reads
hisat2 -p "$THREADS" --dta -x "$REFERENCE" -1 "$r1_fastq" -2 "$r2_fastq" \
	| samtools sort -n -@ "$THREADS" - \
	| samtools fixmate -@ "$THREADS" -m - - \
	| samtools sort -@ "$THREADS" - \
	| samtools markdup -@ "$THREADS" - "$bam_file"
