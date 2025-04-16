#!/bin/bash

set -euo pipefail

BASE_DIR=~/data
RNASEQ_DIR="$BASE_DIR/rnaseq"
TAPESTRI_DIR="$BASE_DIR/tapestri"
DUPES_DIR="$BASE_DIR/duplicates"

mkdir -p "$DUPES_DIR"

# Handle rnaseq -> tapestri
for path in "$RNASEQ_DIR"/*; do
    [[ -d "$path" ]] || continue
    sample=$(basename "$path")

    if [[ "$sample" == *_Mission_Bio ]]; then
        if [[ -d "$TAPESTRI_DIR/$sample" ]]; then
            echo "Duplicate: $sample exists in both. Moving rnaseq copy to duplicates."
            echo mv "$path" "$DUPES_DIR/"
        fi
    fi
done

# Handle tapestri -> rnaseq
for path in "$TAPESTRI_DIR"/*; do
    [[ -d "$path" ]] || continue
    sample=$(basename "$path")

    if [[ "$sample" != *_Mission_Bio ]]; then
        if [[ -d "$RNASEQ_DIR/$sample" ]]; then
            echo "Duplicate: $sample exists in both. Moving tapestri copy to duplicates."
            echo mv "$path" "$DUPES_DIR/"
        fi
    fi
done
