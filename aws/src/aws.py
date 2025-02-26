#!/usr/bin/env python3

import boto3
from botocore.exceptions import ClientError
from util import parse_s3_uri


def list_files(bucket, prefix):
    """List all files in the specified bucket and prefix using boto3."""
    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        operation_params = {"Bucket": bucket, "Prefix": prefix}
        file_list = []
        for page in paginator.paginate(**operation_params):
            if "Contents" in page:
                for obj in page["Contents"]:
                    file_list.append(obj["Key"])
        return file_list
    except ClientError as e:
        print(f"Error listing files: {e}")
        return []


s3_client = boto3.client("s3")

uri1 = "s3://labaf-missionbio-tapestri/AML_singlecell_atlas_JB/"
uri2 = "s3://lab-aaf-ngs-data-archive/RNAseq/10X_Chromium/AML_singlecell_atlas_JB/"

bucket1, prefix1 = parse_s3_uri(uri1)
bucket2, prefix2 = parse_s3_uri(uri2)
# list all files in each dir

files1 = list_files(bucket1, prefix1)
files2 = list_files(bucket2, prefix2,)

# save only the files that end with *.fastq.gz or .md5
fastqs1 = [ file for file in files1 if file.endswith(".fastq.gz") ]
fastqs2 = [ file for file in files2 if file.endswith(".fastq.gz") ]

md5files1 = [ file for file in files1 if file.endswith(".md5") ]
md5files2 = [ file for file in files2 if file.endswith(".md5") ]

# print the number of files in each list
print(f"Found {len(fastqs1)} fastq files in {uri1}")
print(f"Found {len(fastqs2)} fastq files in {uri2}")
print(f"Found {len(md5files1)} md5 files in {uri1}")
print(f"Found {len(md5files2)} md5 files in {uri2}")

# count the total numebr of files and calculate the total size of the files from fastqs1
total_size = 0
for file in fastqs1:
    metadata = s3_client.head_object(Bucket=bucket1, Key=file)
    total_size += metadata["ContentLength"]

print(f"Total size of fastq files in {uri1}: {total_size} bytes")
# convert to TB
print(f"Total size of fastq files in {uri1}: {total_size / 1e12} TB")
