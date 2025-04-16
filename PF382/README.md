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

### Genome Reference

Reference HISAT2 index acquired from `seqs_for_alignment_pipelines.ucsc_ids` at
[RefSeq: NCBI Reference Sequence Database](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/)

## `featureCounts`

Count reads per gene using `featureCounts` from the Subread package.

Output saved as `.tsv` for import into R for DESeq2 analysis.

### Annotation Reference

Reference gtf acquired from `seqs_for_alignment_pipelines.ucsc_ids` at
[RefSeq: NCBI Reference Sequence Database](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/)

## `DESeq2`

Design file, a csv (or any kind of delimiter) file containing:

- filename
- experimental condition
- sample alias

```csv
PF382-CTRL-CLONE-1,CTRL,CTRL_01
PF382-CTRL-CLONE-2,CTRL,CTRL_02
PF382-PHF6-KO-CLONE-2,PHF6,PHF6_02
PF382-PHF6-KO-CLONE-6,PHF6,PHF6_06
PF382-PHF6-KO-CLONE-7,PHF6,PHF6_07
PF382-PHIP-KO-CLONE-19,PHIP,PHIP_19
PF382-PHIP-KO-CLONE-31,PHIP,PHIP_31
PF382-PHIP-KO-CLONE-32,PHIP,PHIP_32
PF382CTRL-CLONE-3,CTRL,CTRL_03
```
