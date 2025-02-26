# RNAseq

Steps:

1. Acquire raw data (FASTQ)
2. Align reads to the reference genome (`hisat2`)
3. Count Reads

## IMPORTANT! PREPARE

Prepare the following qc metrics
- alignment rate
- PCA from deseq2
- quantile normalization of TPM/FPKM

- Try stringtie
  - https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

## FASTQ

S3 URI: `lab-aaf-ngs-data-archive/RNAseq/20241029_Jurkat_PHF6-PHIP-KO_LQ`

```console
JURKAT-CTRL-Clone1-LQ_R1_001.fastq.gz     JURKAT-CTRL-Clone1-LQ_R2_001.fastq.gz
JURKAT-CTRL-Clone2-LQ_R1_001.fastq.gz     JURKAT-CTRL-Clone2-LQ_R2_001.fastq.gz
JURKAT-CTRL-Clone3-LQ_R1_001.fastq.gz     JURKAT-CTRL-Clone3-LQ_R2_001.fastq.gz

JURKAT-PHF6-KO-Clone3-LQ_R1_001.fastq.gz  JURKAT-PHF6-KO-Clone3-LQ_R2_001.fastq.gz
JURKAT-PHF6-KO-Clone7-LQ_R1_001.fastq.gz  JURKAT-PHF6-KO-Clone7-LQ_R2_001.fastq.gz
JURKAT-PHF6-KO-Clone8-LQ_R1_001.fastq.gz  JURKAT-PHF6-KO-Clone8-LQ_R2_001.fastq.gz

JURKAT-PHIP-KO-Clone10-MB_R1_001.fastq.gz JURKAT-PHIP-KO-Clone10-MB_R2_001.fastq.gz
JURKAT-PHIP-KO-Clone22-LQ_R1_001.fastq.gz JURKAT-PHIP-KO-Clone22-LQ_R2_001.fastq.gz
JURKAT-PHIP-KO-Clone5-LQ_R1_001.fastq.gz  JURKAT-PHIP-KO-Clone5-LQ_R2_001.fastq.gz
JURKAT-PHIP-KO-Clone9-MB_R1_001.fastq.gz  JURKAT-PHIP-KO-Clone9-MB_R2_001.fastq.gz
```

## Alignment

GRCh38 (latest major release)
[FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)

Use pre-built index for `hisat2` from `seqs_for_alignment_pipelines.ucsc_ids/`

```sh
# Align reads and process BAM file
# the `dta` flag is only if we're aligning RNAseq reads
hisat2 -p "$THREADS" --dta -x "$REFERENCE" -1 "$r1_fastq" -2 "$r2_fastq" \
        | samtools sort -n -@ "$THREADS" - \
        | samtools fixmate -@ "$THREADS" -m - - \
        | samtools sort -@ "$THREADS" - \
        | samtools markdup -@ "$THREADS" - "$bam_file"
```


## Count Reads

FeatureCounts form Subread package

```sh
featureCounts -a ${ANNOTATION} -p --countReadPairs -T$(nproc) -o total_counts.txt ./*.bam
```

Annotation file from `seqs_for_alignment_pipelines.ucsc_ids/`
- GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz              2024-09-10 13:15   74M  
- GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
