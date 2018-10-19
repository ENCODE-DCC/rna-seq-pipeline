# REFERENCE

This document contains more detailed information on the inputs, outputs and the software.

# CONTENTS

[Software](reference.md#software)  
[Inputs](reference.md#inputs)  
[Outputs](reference.md#outputs)

## Software

### Ubuntu 16.04

The pipeline docker image is based on [Ubuntu base image](https://hub.docker.com/_/ubuntu/) version `16.04`.

### Python 3.5.2

Python parts of the pipeline are run using [Python 3.5.2](https://www.python.org/download/releases/3.5.2/) that ships with Ubuntu 16.04.

### R 3.2.3

MAD QC is run using R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree".

### STAR 2.5.1b

Alignment is done using [STAR 2.5.1b](https://github.com/alexdobin/STAR/releases/tag/2.5.1b). For detailed description of the software see [Article by Dobin et al](https://www.ncbi.nlm.nih.gov/pubmed/23104886). Multiple versions of the software have been released since writing the article.

### RSEM 1.2.23

Quantification is done using [RSEM v1.2.31](https://github.com/deweylab/RSEM/releases/tag/v1.2.31). For detailed description of the software see [Article by Bo Li and Colin N Dewey](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323). Multiple versions of the software have been released since writing the article.

### Kallisto 0.44.0

Additional quantification is calculated using [Kallisto v0.44.0](https://github.com/pachterlab/kallisto/releases/tag/v0.44.0). For detailed description of the software see [Article by Bray et al.](https://www.nature.com/articles/nbt.3519). 

### Samtools 1.9

Samtools flagstats are calculated using [Samtools 1.9](https://github.com/samtools/samtools/releases/tag/1.9). For original publication describing the software, see [Article by Li H et al](https://www.ncbi.nlm.nih.gov/pubmed/19505943). Multiple versions of the software have been released since writing the article.

### bedGraphToBigWig and bedSort

[bedGraphToBigWig kent source version 371](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig) and [bedSort](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort) (no version information available) are used in signal track generation. Source code is downloadable [here](http://hgdownload.soe.ucsc.edu/admin/exe/userApps.v371.src.tgz).
 
## Inputs

A typical input file looks like this:

```
{
    "rna.endedness" : "paired",
    "rna.fastqs_R1" : ["test_data/ENCSR653DFZ_rep1_chr19_10000reads_R1.fastq.gz", "test_data/ENCSR653DFZ_rep2_chr19_10000reads_R1.fastq.gz"],
    "rna.fastqs_R2" : ["test_data/ENCSR653DFZ_rep1_chr19_10000reads_R2.fastq.gz", "test_data/ENCSR653DFZ_rep2_chr19_10000reads_R2.fastq.gz"],
    "rna.aligner" : "star",
    "rna.index" : "test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz",
    "rna.rsem_index" : "test_data/GRCh38_v24_ERCC_phiX_rsemIndex_chr19only.tgz",
    "rna.kallisto.kallisto_index" : "test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx",
    "rna.bamroot" : "PE_stranded",
    "rna.strandedness" : "stranded",
    "rna.strandedness_direction" : "reverse",
    "rna.chrom_sizes" : "test_data/GRCh38_EBV.chrom.sizes",
    "rna.align_ncpus" : 2,
    "rna.align_ramGB" : 4,
    "rna.disks" : "local-disk 20 HDD",
    "rna.kallisto.number_of_threads" : 2,
    "rna.kallisto.ramGB" : 4,
    "rna.rna_qc.tr_id_to_gene_type_tsv" : "transcript_id_to_gene_type_mappings/gencodeV24pri-tRNAs-ERCC-phiX.transcript_id_to_genes.tsv"
}
```

Following elaborates the meaning of each line in the input file.

* `rna.endedness` Indicates whether the endedness of the experiment is `paired` or `single`. 
* `rna.fastqs_R1` Is list of gzipped fastq files containing the first pairs of reads.
* `rna.fastqs_R2` Is list of gzipped fastq files containing the second pairs of reads. 

#### Example:

Assume you are running a paired end experiment with 3 replicates. The fastq files from the first replicate are `replicate1_read1.fastq.gz` and `replicate1_read2.fastq.gz`. The fastq files from the second replicate are `replicate2_read1.fastq.gz` and `replicate2_read2.fastq.gz`. Finally assume that the fastq files from the third replicate are `replicate3_read1.fastq.gz` and `replicate3_read2.fastq.gz`. In this case the input on the relevant part should be as follows:  
`"rna.fastqs_R1" : ["replicate1_read1.fastq.gz", "replicate2_read1.fastq.gz", "replicate3_read1.fastq.gz"]`  
`"rna.fastqs_R2" : ["replicate1_read2.fastq.gz", "replicate2_read2.fastq.gz", "replicate3_read2.fastq.gz"]`  
Note that it is very important that the replicates are in same order in both lists, this correspondence is used for pairing correct files with each other.

* `rna.aligner` Use `star` aligner, possibly extended to use others in future.
* `rna.index` Is the index for STAR aligner. 
* `rna.rsem_index` Is the index for RSEM quantifier.
* `rna.kallisto.kallisto_index` Is the index for Kallisto quantifier.
* `rna.bamroot` This is a prefix that gets added into the output filenames. Additionally the files are prefixed with information of the replicate they originate from.

#### Example:

Assume the `rna.bamroot` is `FOO`. Outputs from first replicate would be prefixed by `rep1FOO` and outputs from second replicate would be prefixed by `rep2FOO` etc.

* `rna.strandedness` Indicates whether the experiment is `stranded` or `unstranded`.
* `rna.strandedness_direction` Indicates the direction of strandedness. Options are `forward`, `reverse` and `unstranded`.
* `rna.chrom_sizes` Is the file containing the chromosome sizes. You can find and download the files from [ENCODE portal](https://www.encodeproject.org/references/ENCSR425FOI/).
* `rna.align_ncpus` How many cpus are available for STAR alignment and RSEM quantification.
* `rna.align_ramGB` How many GBs of memory are available for STAR alignment and RSEM quantification.
* `rna.disks` How much disk space is available for pipeline. You can also specify the type of disk, `HDD` for a spinning disk and `SSD` for a solid state drive.
* `rna.kallisto.number_of_threads` How many cpus are available for Kallisto quantification.
* `rna.kallisto.ramGB` How many GBs of memory are available for Kallisto quantification.

#### Example:

Assume you want to allocate 100 gigabytes of spinning hard drive. In this case you would enter `"local-disk 100 HDD"`. If you want to allocate 111 gigabytes of solid state drive space, enter `"local-disk 111 SSD"`.

#### Additional inputs when running single-ended experiments:

Kallisto quantifier makes use of average fragment lenghts and standard deviations of those lengths. In the case of paired end experiments, those values can be calculated from the data, but in case of single-ended experiment those values must be provided.

* `rna.kallisto.fragment_length` Is the average fragment length.
* `rna.kallisto.sd_of_fragment_length` Is the standard deviation of the fragment lengths.
* `rna.rna_qc.tr_id_to_gene_type_tsv` rna_qc task calculates the number of reads by gene type. For this a tsv file that contains a mapping from transcript IDs to gene types is provided. For GRCh38, hg19, and mm10 with ERCC (ambion 1) and PhiX spikes the tsv is provided in this repo. If you are using some other annotation, you can use code [here](https://github.com/ENCODE-DCC/transcript_id_to_gene_type_mapping) to build your own.

## Outputs

1. `DNANexus`: If you choose to use `dxWDL` and run pipelines on DNANexus platform, then output will be stored on the specified output directory without any subdirectories.

2. `Cromwell`: `Cromwell` will store outputs for each task under directory `cromwell-executions/[WORKFLOW_ID]/call-[TASK_NAME]/shard-[IDX]`. For all tasks `[IDX]` means a zero-based index for each replicate.

### Output files

#### Task Align

* `Genome bam`, file name ends with `_genome.bam`. Bam aligned to genome.
* `Anno bam`, file name ends with `_anno.bam`. Bam aligned to annotation.
* `Genome flagstat` file name ends with `_genome_flagstat.txt`. Samtools flagstats on the genome bam.
* `Anno flagstat` file name ends with `_anno_flagstat.txt`. Samtools flagstats on anno bam.
* `STAR run log` file name ends with `_Log.final.out`. STAR run log.
* `Python log` file name is `align.log`. This file contains possible additional information on the pipeline step.

#### Task Kallisto

* `Kallisto quants`, file name ends with `_abundance.tsv`. Kallisto quantifications.
* `Python log` file name is `kallisto_quant.log`. This file contains possible additional information on the pipeline step.

#### Task Bam to Signals

In case of an stranded run, the plus and minus strand signal tracks are separated (there will be four tracks per replicate).

* `Unique BigWig`, file name ends with `niq.bw`. Contains the signal track of the uniquely mapped reads.
* `All BigWig`, the file name ends with `ll.bw`. Contains the signal track of all reads.
* `Python log` file name is `bam_to_signals.log`. This file contains possible additional information on the pipeline step.

#### Task RSEM Quant

* `Genes results`, file name ends with `genes.results`. Contains gene quantifications.
* `Isoforms results`, file name ends with `isoforms.results`. Contains isoform quantifications.
* `Number of genes`, file name ends with `_number_of_genes_detected.json`. Contains the number of genes detected, which is determined as `TPM` value being greater than `1`.
* `Python log` file name is `rsem_quant.log`. This file contains possible additional information on the pipeline step.

### Task Mad QC

This step is run if and only if the number of replicates is 2.

* `Mad QC plot`, file name ends with `_mad_plot.png`. Contains the MAD QC plot.
* `Mad QC metrics` file name ends with `_mad_qc_metrics.json`. Contains MAD QC metrics.
* `Python log` file name is `mad_qc.log`. This file contains possible additional information on the pipeline step.

### Task RNA QC

This step calculates additional metrics. At this time the only metric is to calculate reads by gene type. It is very **IMPORTANT** to look at the `Python log` of this step to see that the transcriptome bam did not contain any transcripts that are not present in the transcript ID to gene type mapping tsv. In case that happens, make sure you are using the STAR aligner and RSEM quantifier indexes you think you are using, and that all the other references are correct!

* `RNA QC`, file name ends with `qc.json`. Contains additional QC metrics. For now the reads by gene type.
* `Python log` file name is `rna_qc.log`. This file contains **IMPORTANT** information on the pipeline step.
