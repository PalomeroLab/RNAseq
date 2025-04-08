#!/bin/bash
# Set up some environment variables for the PF382 project
ASSAY=RNAseq
EXPERIMENT=20250310_PF382_PHF6-PHIP-KO_LQ
export S3URI="s3://lab-aaf-ngs-data-archive/${ASSAY}/"${EXPERIMENT}
