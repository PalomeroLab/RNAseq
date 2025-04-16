#!/usr/bin/env bash

set -eu

LOCALDIR=/mnt/data
S3URI="s3://labaf-missionbio-tapestri"

while IFS= read -r line; do
	echo "Syncing $line"
	sudo aws s3 cp "$S3URI/$line" "$LOCALDIR/$line" --force-glacier-transfer
	# done < ./md5files1.txt
done < ./fastqs1.txt

sudo aws s3 sync s3://labaf-missionbio-tapestri/AML_singlecell_atlas_JB/ /mnt/data/AML_singlecell_atlas_JB/ --size-only --exclude "*" --include "*.fastq.gz" --force-glacier-transfer 
