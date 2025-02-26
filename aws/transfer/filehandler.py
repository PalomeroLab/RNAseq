#!/usr/bin/env python3

import polars as pl

class FastqFile:
    def __init__(self, filename):
        self.filename = filename
        self.md5 = False

class Sample:
    def __init__(self, dirname, sampleID, fastqfiles):
        self.dirname = dirname
        self.sampleID = sampleID
        self.fastqfiles = fastqfiles

def process_files(fastq_file_path, md5_file_path):
    with open(fastq_file_path, "r") as f:
        fastq_files = [line.strip() for line in f.readlines()]

    with open(md5_file_path, "r") as f:
        md5_files = [line.strip() for line in f.readlines()]

    samples = []
    seen_samples = set()  # Track seen (dirname, sampleID) combinations

    for fastq_file_path in fastq_files:
        parts = fastq_file_path.split("/")
        if len(parts) < 3:
            continue

        dir_name = parts[0]
        sample_id = parts[1]
        sample_file_name = parts[-1]  # Corrected filename extraction

        key = (dir_name, sample_id)
        if key not in seen_samples:
            new_sample = Sample(dir_name, sample_id, [])
            samples.append(new_sample)
            seen_samples.add(key)

        fastq_obj = FastqFile(sample_file_name) # Pass only filename
        # Find the Sample object and append (efficient because of seen_samples)
        for sample in samples:
            if (sample.dirname, sample.sampleID) == key:
                sample.fastqfiles.append(fastq_obj)
                break  # Important: Stop searching once found


    for md5_file_path in md5_files:
        md5_filename = md5_file_path.split("/")[-1]
        for sample in samples:
            for fastq_obj in sample.fastqfiles:
                if md5_filename.replace(".md5", "") in fastq_obj.filename:
                    fastq_obj.md5 = True
                    break
            else:
                continue
            break

    return samples

def main():
    fastq_file_path = "sync/fastqs1.txt"
    md5_file_path = "sync/md5files1.txt"
    samples = process_files(fastq_file_path, md5_file_path)

    for sample in samples:
        print(f"Directory: {sample.dirname}")
        print(f"Sample ID: {sample.sampleID}")
        print("Fastq Files:")
        for fq_file_obj in sample.fastqfiles:
            print(f"  - {fq_file_obj.filename} (MD5 Exists: {fq_file_obj.md5})")
        print("-" * 20)

    data_for_polars = []
    for sample in samples:
        for fastq_obj in sample.fastqfiles:
            data_for_polars.append({
                "directory": sample.dirname,
                "sample": sample.sampleID,
                "filename": fastq_obj.filename,
                "md5_exists": fastq_obj.md5,
            })

    df = pl.DataFrame(data_for_polars)
    output_csv_path = "output.csv"
    df.write_csv(output_csv_path)
    print(f"DataFrame written to {output_csv_path}")

if __name__ == "__main__":
    main()
