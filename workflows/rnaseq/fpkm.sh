#!/bin/sh

CTS=./data/LAURA/PHIP/counts_noclone9.txt

./fpkmFromFeatureCounts.R ${1:-$CTS}

