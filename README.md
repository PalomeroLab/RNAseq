# RNA Sequencing Analysis

RNA sequencing (RNA-seq) is a technique used to analyze the _transcriptome_ —
the complete set of RNA transcripts in a cell. This method involves sequencing
RNA molecules after reverse transcription to cDNA to examine gene expression levels
and identify novel transcripts.

## Pre-requisites

Before starting the analysis, ensure you have the following:

- Aligned reads in SAM/BAM/CRAM format
- A reference annotation file in GTF or GFF format
- A read quantification tool such as featureCounts or HTSeq

This document describes two pipelines for RNA-seq analysis:

- Using `featureCounts` and then analyzing with DESeq2 in R
- Using the Tuxedo Suite (HISAT2, StringTie, Ballgown)

## `featureCounts`

[`featureCounts`](https://subread.sourceforge.net/featureCounts.html)
is a program for counting reads mapped to genomic features, such as genes, exons,
and promoters. It is part of the [Subread](https://subread.sourceforge.net)
package and can be used for RNA-seq as well as DNA-seq analysis.

> featureCounts takes as input SAM/BAM files and an annotation file including
> chromosomal coordinates of features. It outputs numbers of reads assigned to features
> (or meta-features). It also outputs stat info for the overall summrization results,
> including number of successfully assigned reads and number of reads that failed to be
> assigned due to various reasons (these reasons are included in the stat info).

Use FeatureCounts to count the number of reads that map to each gene in a GTF
file and summarize the results for downstream analysis (i.e., differential expression).

### Input

Aligned reads in BAM format and a GTF file containing genomic features.

Example usage:

```sh
featureCounts -T "$NUM_THREADS" --verbose -t exon -g gene_id --countReadPairs \
  -a "$REF_GTF" -p -P -C -B -o "${OUT_DIR}/${OUT_PREFIX}.tsv" ./*.bam
```

Or, more simply:

```sh
featureCounts -a ${ANNOTATION} -p --countReadPairs -T$(nproc) -o total_counts.txt ./*.bam
```

Annotation file from `seqs_for_alignment_pipelines.ucsc_ids/`

- GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz 2024-09-10 13:15 74M
- GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz

### Output

Feature counts are written to a tab-delimited file (`.tsv` or `.txt`)
with columns for each sample. You can import this file into R or other
statistical software for further analysis.

Feature counts also provide a summary of the number of reads that were
assigned to features, as well as the number of reads that were not assigned
to any feature.

### Downstream Analysis using R

Export the results from `featureCounts` to R for further analysis, such as
differential expression analysis using packages like DESeq2 or edgeR.

## Normalization

> [!CAUTION]
> Without between sample normalization NONE of these units are comparable across experiments.
> This is a result of RNA-Seq being a relative measurement, not an absolute one.
> [Source](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).

**Counts** usually refers to the number of reads that align to a particular feature ($X_i).
These numbers are heavily dependent on two things:

- The amount of fragments you sequenced (this is related to relative abundances).
- The length of the feature, or more appropriately, the _effective length_.

Effective length refers to the number of possible start sites a feature could
have generated a fragment of that particular length.

Since counts are NOT scaled by the length of the feature, all units in this
category are not comparable within a sample without adjusting for the feature
length. This means you can’t sum the counts over a set of features to get the
expression of that set (e.g. you can’t sum isoform counts to get gene counts).

### TPM (Transcripts Per Million)

Transcripts per million (TPM) is a measurement of the proportion of transcripts
in your pool of RNA.

### RPKM/FPKM

Reads per kilobase of exon per million reads mapped (RPKM), or the more generic
FPKM (substitute reads with fragments) are essentially the same thing.

> [!NOTE]
> FPKM is not $2 *$ RPKM if you have paired-end reads. FPKM $==$ RPKM if you have single-end reads,
> and saying RPKM when you have paired-end reads is incorrect.

<!-- TODO: add formulae -->

### Code example

```R
countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}
```

## R Packages

### DESeq2

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is an
R package for differential gene expression analysis based on the negative
binomial distribution.

> DESeq2 provides methods to test for differential expression by use of negative
> binomial generalized linear models. The models use the raw counts as input and
> perform regularized log transformation and variance stabilizing
> transformation. The package also provides functions to visualize the data and
> results.

## Tuxedo Suite

The Tuxedo Suite is a collection of tools for transcript-level expression analysis of RNA-seq experiments.

1. Align reads to the reference genome and sort to BAM format (HISAT2)
2. Assemble transcripts (StringTie)
3. Prepare for differential expression analysis (Ballgown setup)
4. Perform differential expression analysis
5. Visualize results

| Tool      | Description                                    | Manual                                                                      | Source                                             |
| --------- | ---------------------------------------------- | --------------------------------------------------------------------------- | -------------------------------------------------- |
| HISAT2    | Align reads to reference genome                | [manual](https://daehwankimlab.github.io/hisat2/manual/)                    | [source](https://github.com/DaehwanKimLab/hisat2)  |
| StringTie | Assemble RNA-Seq alignments into transcripts   | [manual](https://ccb.jhu.edu/software/stringtie/index.shtml)                | [source](https://github.com/gpertea/stringtie)     |
| Ballgown  | Isoform-level differential expression analysis | [manual](https://bioconductor.org/packages/release/bioc/html/ballgown.html) | [source](https://github.com/alyssafrazee/ballgown) |

> [!TIP]
> StringTie comes packaged with `gffcompare` for comparing and evaluating the
> accuracy of RNA-seq transcript assemblers. Read the
> [manual](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) for details.

Using Ballgown:

1. Filter to remove low-abundance genes
2. Identify differentially expressed transcripts and genes
3. Add gene names
4. Sort by p-value
5. Write results to file

## References

- Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-30. [PMID: 24227677](https://pubmed.ncbi.nlm.nih.gov/24227677/)
- Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. [PMID: 27560171](https://pubmed.ncbi.nlm.nih.gov/27560171/)
- Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357-360. [PMID: 25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)
- Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015;33(3):290-295. [PMID: 25690850](https://pubmed.ncbi.nlm.nih.gov/25690850/)
- Frazee AC, Pertea G, Jaffe AE, Langmead B, Salzberg SL, Leek JT. Ballgown bridges the gap between transcriptome assembly and expression analysis. Nat Biotechnol. 2015;33(3):243-246. [PMID: 25748911](https://pubmed.ncbi.nlm.nih.gov/25748911/)
