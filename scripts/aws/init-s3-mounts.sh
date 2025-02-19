#!/bin/bash
## Mounts S3 buckets to local directories
## https://www.petergirnus.com/blog/how-to-mount-an-aws-s3-bucket-on-macos
## Requires:
##	s3fs (on macos)
##	s3-mount (on linux)

S3_MOUNTDIR="${HOME}/S3"

# Stick to bash for some array support
S3_BUCKETS=(
	kalay
	lab-aaf-ngs-data-archive
	lab-aaf-scratch
	lab-aaf-server-backup
	lab-tp-inventory
	lab-tp-rstudio-scratch
	labaf-missionbio-tapestri
	lq2230
	palomerolab
	ptcl-ctcl-gene-expr-database
)

for bucket in "${S3_BUCKETS[@]}"; do
	mkdir -p "${S3_MOUNTDIR}/${bucket}"
	# s3-mount --allow-other "${bucket}" "${S3_MOUNTDIR}/${bucket}"
	# NOTE: use s3fs instead of s3-mount if on macos
	s3fs "${bucket}" "${S3_MOUNTDIR}/${bucket}"
	echo "Mounted ${bucket} to ${S3_MOUNTDIR}/${bucket}"
done
