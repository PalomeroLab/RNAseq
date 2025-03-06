#!/bin/sh
set -eux

ANNOTATION="/opt/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz"

FLAGS="-a ${ANNOTATION} -p --countReadPairs -T$(nproc)"

featureCounts "$FLAGS" -o counts.txt ./*.bam
