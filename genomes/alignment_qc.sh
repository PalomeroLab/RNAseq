#!/bin/sh

# Run all the samtools commands supported by multiqc:
# - stats
# - flagstats
# - idxstats
# - rmdup
# - coverage
# - markdup
set -eux

# First param is the folder containing the .bam files
if [ $# -lt 1 ]; then
	echo "Usage: $0 <bam_folder>"
	exit 1
fi

OUTDIR="$1"/../samtools_output
mkdir -p "$OUTDIR"

for f in ./"$1"/*.bam; do
	# samtools stats -@"$(nproc)" "$f" > "$OUTDIR"/"$(basename "$f" .bam)".stats
	# samtools flagstat -@"$(nproc)" "$f" > "$OUTDIR"/"$(basename "$f" .bam)".flagstat
	samtools coverage  "$f" > "$OUTDIR"/"$(basename "$f" .bam)".coverage &
done
