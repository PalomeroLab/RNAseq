#!/usr/bin/env python3

import os

d1 = "/Users/rdn/GitHub/palomerolab/aws/aml/md5s1"
d2 = "/Users/rdn/GitHub/palomerolab/aws/aml/md5s2"

# go throuh each file in d1 and compare to d2
# not any files that are not identical
# not any files that are in d2 but not in d1 and vice versa

files1 = os.listdir(d1)
files2 = os.listdir(d2)

for file in files1:
    if file not in files2:
        print(f"{file} not in {d2}")

for file in files2:
    if file not in files1:
        print(f"{file} not in {d1}")

for file in files1:
    if file in files2:
        with open(f"{d1}/{file}") as f1:
            with open(f"{d2}/{file}") as f2:
                if f1.read() != f2.read():
                    print(f"{file} not identical")


# count the number of files in each directory
print(f"Found {len(files1)} files in {d1}")
print(f"Found {len(files2)} files in {d2}")

# compares the files in
# check files in ./fastqs{1,2}.txt and note any not in both
fastqs1 = []
fastqs2 = []
with open("fastqs1.txt") as f:
    fastqs1 = f.readlines()
with open("fastqs2.txt") as f:
    fastqs2 = f.readlines()

# compare the two lists: 
# RNAseq/10X_Chromium/AML_singlecell_atlas_JB/CU_M_115_Mission_Bio/JM75_R1_001.fastq.gz
# AML_singlecell_atlas_JB/CU_M_185_Dx/fastqs/AJ027_S1_L001_R1_001.fastq.gz
# remove `RNAseq/10X_Chromium/AML_singlecell_atlas_JB/` or  `AML_singlecell_atlas_JB`
fastqs1 = [ file.split("/")[-1].strip() for file in fastqs1 ]
fastqs2 = [ file.split("/")[-1].strip() for file in fastqs2 ]

for file in fastqs1:
    if file not in fastqs2:
        print(f"{file} not in fastqs2.txt")

for file in fastqs2:
    if file not in fastqs1:
        print(f"{file} not in fastqs1.txt")
