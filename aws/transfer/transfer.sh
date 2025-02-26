#!/bin/bash
set -eux

SOURCEROLE=
# Create an IAM user and role in the destination AWS account
aws iam create-role \
	--role-name S3MigrationRole \
	--assume-role-policy-document file://./S3MigrationRole.json

# Create and attach the S3 bucket policy in the source account
aws s3api put-bucket-policy --bucket labaf-missionbio-tapestri --policy file://./S3BucketPolicy.json


aws sts assume-role --role-arn "arn:aws:iam::618638526287:role/S3MigrationRole" --role-session-name AWSCLI-Session

aws sts assume-role \
    --role-arn "arn:aws:iam::925204023671:role/S3MigrationRole" \
    --role-session-name AWSCLI-Session


aws s3 cp s3://lab-aaf-ngs-data-archive/RNAseq/10X_Chromium/AML_singlecell_atlas_JB/README.txt \
    s3://ziwei-collaborator-upload --profile aml

aws s3 cp s3://lab-aaf-ngs-data-archive/RNAseq/10X_Chromium/AML_singlecell_atlas_JB/README.txt \
