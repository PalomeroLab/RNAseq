#!/bin/bash
## run all the samtools commands supported by multiqc:
## - stats
## - flagstats
## - idxstats
## - rmdup
## - coverage
## - markdup
set -eux

# first param is the folder containing the .bam files
if [ $# -lt 1 ]; then
	echo "usage: $0 <bam_folder>"
	exit 1
fi

outdir="$1"/../samtools_output
mkdir -p "$outdir"

for f in ./"$1"/*.bam; do
	# samtools stats -@"$(nproc)" "$f" > "$outdir"/"$(basename "$f" .bam)".stats
	# samtools flagstat -@"$(nproc)" "$f" > "$outdir"/"$(basename "$f" .bam)".flagstat
	samtools coverage  "$f" > "$outdir"/"$(basename "$f" .bam)".coverage &
done
