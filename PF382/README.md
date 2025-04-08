# RNAseq/20250310_PF382_PHF6-PHIP-KO_LQ/

Sequenced by Azenta

## Manifest

```sh
aws s3 ls s3://lab-aaf-ngs-data-archive/RNAseq/20250310_PF382_PHF6-PHIP-KO_LQ/00_fastq/
```

- PF382-CTRL-CLONE-1_R1_001.fastq.gz
- PF382-CTRL-CLONE-1_R2_001.fastq.gz
- PF382-CTRL-CLONE-2_R1_001.fastq.gz
- PF382-CTRL-CLONE-2_R2_001.fastq.gz
- PF382-PHF6-KO-CLONE-2_R1_001.fastq.gz
- PF382-PHF6-KO-CLONE-2_R2_001.fastq.gz
- PF382-PHF6-KO-CLONE-6_R1_001.fastq.gz
- PF382-PHF6-KO-CLONE-6_R2_001.fastq.gz
- PF382-PHF6-KO-CLONE-7_R1_001.fastq.gz
- PF382-PHF6-KO-CLONE-7_R2_001.fastq.gz
- PF382-PHIP-KO-CLONE-19_R1_001.fastq.gz
- PF382-PHIP-KO-CLONE-19_R2_001.fastq.gz
- PF382-PHIP-KO-CLONE-31_R1_001.fastq.gz
- PF382-PHIP-KO-CLONE-31_R2_001.fastq.gz
- PF382-PHIP-KO-CLONE-32_R1_001.fastq.gz
- PF382-PHIP-KO-CLONE-32_R2_001.fastq.gz
- PF382CTRL-CLONE-3_R1_001.fastq.gz
- PF382CTRL-CLONE-3_R2_001.fastq.gz

> [!NOTE]
> Typo in the last sample name; should fix it in the future.

## Alignment

Align paired end samples against the human genome (hg38) using HISAT2.

Pipeline to `.bam`:

- `fixmated`'d
- sorted
- duplicates marked (but not removed)

## `featureCounts`

Count reads per gene using `featureCounts` from the Subread package.

Output saved as `.tsv` for import into R for DESeq2 analysis.

## `DESeq2`

WIP...
