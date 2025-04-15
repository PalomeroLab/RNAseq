#!/bin/bash

CURRENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load the environment
. "$CURRENT_DIR"/../setup.sh

aws s3 sync "$S3URI"/02_counts/ .

# md5sum -c md5sums.txt
