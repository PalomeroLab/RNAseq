#!/usr/bin/env python3

import os
import sys
import polars as pl


class FastqFile:
    def __init__(self, filename):
        self.filename = filename
        self.md5 = False


class Sample:
    def __init__(self, dirname, sample_id):
        self.dirname = dirname
        self.sample_id = sample_id
        self.fastqfiles = []


def collect_files(directory):
    fastq_paths = []
    md5_paths = []
    for root, _, files in os.walk(directory):
        for f in files:
            full_path = os.path.join(root, f)
            if f.endswith(".fastq.gz"):
                fastq_paths.append(full_path)
            elif f.endswith(".md5"):
                md5_paths.append(full_path)
    return fastq_paths, md5_paths


def build_samples(fastq_paths, md5_paths):
    samples = {}
    for fq in fastq_paths:
        parts = fq.strip().split("/")
        if len(parts) < 3:
            continue
        dirname, sample_id = parts[-3], parts[-2]
        filename = parts[-1]
        key = (dirname, sample_id)
        if key not in samples:
            samples[key] = Sample(dirname, sample_id)
        samples[key].fastqfiles.append(FastqFile(filename))

    for md5 in md5_paths:
        md5_name = os.path.basename(md5).replace(".md5", "")
        for sample in samples.values():
            for fq in sample.fastqfiles:
                if md5_name in fq.filename:
                    fq.md5 = True
    return samples.values()


def report_samples(samples):
    records = []
    for s in samples:
        print(f"Directory: {s.dirname}\nSample ID: {s.sample_id}")
        for fq in s.fastqfiles:
            print(f"  - {fq.filename} (MD5 Exists: {fq.md5})")
            records.append(
                {
                    "directory": s.dirname,
                    "sample": s.sample_id,
                    "filename": fq.filename,
                    "md5_exists": fq.md5,
                }
            )
        print("-" * 20)
    return records


def main():
    if len(sys.argv) != 2:
        print("Usage: python filehandler.py <directory>")
        sys.exit(1)
    directory = sys.argv[1]
    fastq_paths, md5_paths = collect_files(directory)
    samples = build_samples(fastq_paths, md5_paths)
    data = report_samples(samples)
    df = pl.DataFrame(data)
    df.write_csv("output.csv")


if __name__ == "__main__":
    main()
