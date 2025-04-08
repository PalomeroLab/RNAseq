#!/bin/bash

CURRENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load the environment
. "${CURRENT_DIR:-.}"/../setup.sh

aws s3 sync . "$S3URI"/reports/ --exclude "*.sh"

