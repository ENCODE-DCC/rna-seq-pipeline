# REFERENCE

This document contains more detailed information on the inputs, outputs and the software.

# CONTENTS

[Software](reference.md#software)  
[Genome Reference Files](reference.md#genome-reference-files)  
[Inputs](reference.md#inputs)  
[Resource Considerations](reference.md#note-about-resources)  
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

### RSEM 1.2.31

Quantification is done using [RSEM v1.2.31](https://github.com/deweylab/RSEM/releases/tag/v1.2.31). For detailed description of the software see [Article by Bo Li and Colin N Dewey](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323). Multiple versions of the software have been released since writing the article.

### Kallisto 0.44.0

Additional quantification is calculated using [Kallisto v0.44.0](https://github.com/pachterlab/kallisto/releases/tag/v0.44.0). For detailed description of the software see [Article by Bray et al.](https://www.nature.com/articles/nbt.3519).

### Samtools 1.9

Samtools flagstats are calculated using [Samtools 1.9](https://github.com/samtools/samtools/releases/tag/1.9). For original publication describing the software, see [Article by Li H et al](https://www.ncbi.nlm.nih.gov/pubmed/19505943). Multiple versions of the software have been released since writing the article.

### bedGraphToBigWig and bedSort

[bedGraphToBigWig kent source version 371](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig) and [bedSort](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort) (no version information available) are used in signal track generation. Source code is downloadable [here](http://hgdownload.soe.ucsc.edu/admin/exe/userApps.v371.src.tgz).

## Genome Reference Files

Reference and index files can be downloaded from the [ENCODE Portal](https://www.encodeproject.org/search/?type=Reference&reference_type=index). These files are for human, and mouse will follow. Files that are needed are STAR Index, RSEM Index, Kallisto Index, chromosome sizes file. Links to required files are summarized in the table below.

|Genome|Annotation Version|STAR|RSEM|Kallisto|Chr Sizes|
|--------|:--------------:|----|----|--------|---------|
|[GRCh38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)|Gencode V29|[GRCh38 V29 STAR index](https://www.encodeproject.org/files/ENCFF598IDH/)|[GRCh38 V29 RSEM index](https://www.encodeproject.org/files/ENCFF285DRD/)|[V29 Kallisto index](https://www.encodeproject.org/files/ENCFF471EAM/)|[GRCh38 Chr sizes](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)|
|[GRCh38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)|Gencode V24|[GRCh38 V24 STAR index](https://www.encodeproject.org/files/ENCFF742NER/)|[GRCh38 V24 RSEM index](https://www.encodeproject.org/references/ENCSR219BJA/)|[V24 Kallisto index](https://www.encodeproject.org/references/ENCSR622RMG/)|[GRCh38 Chr sizes](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)|
|[mm10](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/)|Gencode M21|[mm10 M21 STAR index](https://www.encodeproject.org/files/ENCFF795ZFF/)|[mm10 M21 RSEM index](https://www.encodeproject.org/files/ENCFF363TFV/)|[M21 Kallisto index](https://www.encodeproject.org/files/ENCFF383TUX/)|[mm10 Chr sizes](https://www.encodeproject.org/files/mm10_no_alt.chrom.sizes/)|

## Inputs

Input is given in a json file like this (note that the resources in this input file are suitable for small subsampled test files, for recommendations that work for typical full sized files, see [Resource Considerations](reference.md#note-about-resources)):

```
{
    "rna.endedness" : "paired",
    "rna.fastqs_R1" : [["test_data/ENCSR653DFZ_rep1_chr19_10000reads_R1.fastq.gz"], ["test_data/ENCSR653DFZ_rep2_chr19_10000reads_R1.fastq.gz"]],
    "rna.fastqs_R2" : [["test_data/ENCSR653DFZ_rep1_chr19_10000reads_R2.fastq.gz"], ["test_data/ENCSR653DFZ_rep2_chr19_10000reads_R2.fastq.gz"]],
    "rna.aligner" : "star",
    "rna.align_index" : "test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz",
    "rna.rsem_index" : "test_data/GRCh38_v24_ERCC_phiX_rsemIndex_chr19only.tgz",
    "rna.bamroot" : "PE_stranded",
    "rna.strandedness" : "stranded",
    "rna.strandedness_direction" : "reverse",
    "rna.chrom_sizes" : "test_data/GRCh38_EBV.chrom.sizes",
    "rna.align_ncpus" : 2,
    "rna.align_ramGB" : 4,
    "rna.kallisto_index" : "test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx",
    "rna.kallisto_number_of_threads" : 2,
    "rna.kallisto_ramGB" : 4,
    "rna.rna_qc_tr_id_to_gene_type_tsv" : "transcript_id_to_gene_type_mappings/gencodeV24pri-tRNAs-ERCC-phiX.transcript_id_to_genes.tsv",
    "rna.bam_to_signals_ncpus" : 1,
    "rna.bam_to_signals_ramGB" : 2,
    "rna.rsem_ncpus" : 2,
    "rna.rsem_ramGB" : 4,
    "rna.align_disk" : "local-disk 20 HDD",
    "rna.kallisto_disk" : "local-disk 20 HDD",
    "rna.rna_qc_disk" : "local-disk 20 HDD",
    "rna.mad_qc_disk" : "local-disk 20 HDD",
    "rna.bam_to_signals_disk" : "local-disk 20 HDD",
    "rna.rsem_disk" : "local-disk 20 HDD"
}
```

Following elaborates the meaning of each line in the input file.

* `rna.endedness` Indicates whether the endedness of the experiment is `paired` or `single`.
* `rna.fastqs_R1` Is list of lists of gzipped fastq files containing the first pairs of reads.
* `rna.fastqs_R2` Is list of lists of gzipped fastq files containing the second pairs of reads.

#### Example:

Assume you are running a paired end experiment with 3 replicates. The fastq files from the first replicate are `replicate1_read1.fastq.gz` and `replicate1_read2.fastq.gz`. The fastq files from the second replicate are `replicate2_read1.fastq.gz` and `replicate2_read2.fastq.gz`. Finally assume that the fastq files from the third replicate are `replicate3_read1.fastq.gz` and `replicate3_read2.fastq.gz`. In this case the input on the relevant part should be as follows:
`"rna.fastqs_R1" : [["replicate1_read1.fastq.gz"], ["replicate2_read1.fastq.gz"], ["replicate3_read1.fastq.gz"]]`
`"rna.fastqs_R2" : [["replicate1_read2.fastq.gz"], ["replicate2_read2.fastq.gz"], ["replicate3_read2.fastq.gz"]]`
Note that it is very important that the replicates are in same order in both lists, this correspondence is used for pairing correct files with each other.

#### Example:

Assume you are running a single end experiment with 3 replicates. Further assume that in your first round of sequencing you do not get depth you would like (but you still think that the data you obtained is good), and that the data for the replicates are in files `replicate1_part1.fastq.gz` and `replicate2_part1.fastq.gz`. To obtain a deeper data, you submit the libraries from the replicates to re-sequencing, which yields data files  `replicate1_part2.fastq.gz` and `replicate2_part2.fastq.gz`. In this case the input on the fastq part should be:
`"rna.fastqs_R1" : [["replicate1_part1.fastq.gz", "replicate1_part2.fastq.gz"], ["replicate2_part1.fastq.gz", "replicate2_part2.fastq.gz"]]`
`"rna.fastqs_R2" : []` (in single ended experiment the second read input is empty)
In the case when the input for a replicate consists of several fastq files, they will be automatically merged during the pipeline run.

* `rna.aligner` Use `star` aligner, possibly extended to use others in future.
* `rna.align_index` Is the index for STAR aligner.
* `rna.rsem_index` Is the index for RSEM quantifier.
* `rna.kallisto_index` Is the index for Kallisto quantifier.
* `rna.bamroot` This is a prefix that gets added into the output filenames. Additionally the files are prefixed with information of the replicate they originate from.

#### Example:

Assume the `rna.bamroot` is `FOO`. Outputs from first replicate would be prefixed by `rep1FOO` and outputs from second replicate would be prefixed by `rep2FOO` etc.

* `rna.strandedness` Indicates whether the experiment is `stranded` or `unstranded`. If this is `stranded`, then the `rna.strandedness_direction` should be set to `forward` or `reverse`.
* `rna.strandedness_direction` Indicates the direction of strandedness. Options are `forward`, `reverse` and `unstranded`.
* `rna.chrom_sizes` Is the file containing the chromosome sizes. You can find and download the files from ENCODE Portal, see [table above](reference.md#genome-reference-files).
* `rna.align_ncpus` How many cpus are available for STAR alignment.
* `rna.align_ramGB` How many GBs of memory are available for STAR alignment.
* `rna.align_ncpus` How many cpus are available for RSEM quantification.
* `rna.align_ramGB` How many GBs of memory are available for RSEM quantification.
* `rna.align_disk` How much disk space is available for Align task. You can also specify the type of disk, `HDD` for a spinning disk and `SSD` for a solid state drive.
* `rna.kallisto_disk` As above, but for Kallisto.
* `rna.rna_qc_disk` As above, but for RNA QC.
* `rna.bam_to_signals_disk` As above, but for bam_to_signals.
* `rna.mad_qc_disk` As above, but for MAD QC.
* `rna.rsem_disk` As above, but for RSEM.
* `rna.kallisto_number_of_threads` How many threads are available for Kallisto quantification.
* `rna.kallisto_ramGB` How many GBs of memory are available for Kallisto quantification.


#### Example:

Assume you want to allocate 100 gigabytes of spinning hard drive. In this case you would enter `"local-disk 100 HDD"`. If you want to allocate 111 gigabytes of solid state drive space, enter `"local-disk 111 SSD"`.

* `rna.rna_qc_tr_id_to_gene_type_tsv` rna_qc task calculates the number of reads by gene type. For this a tsv file that contains a mapping from transcript IDs to gene types is provided. For GRCh38, hg19, and mm10 with ERCC (ambion 1) and PhiX spikes the tsv is provided in this repo. If you are using some other annotation, you can use code [here](https://github.com/ENCODE-DCC/transcript_id_to_gene_type_mapping) to build your own.
* `rna.bam_to_signals_ncpus` Is the number of cpus given to bam_to_signals task.
* `rna.bam_to_signals_ramGB` Is the amount of memory in GB given to bam_to_signals task.

#### Additional inputs when running single-ended experiments:

Kallisto quantifier makes use of average fragment lenghts and standard deviations of those lengths. In the case of paired end experiments, those values can be calculated from the data, but in case of single-ended experiment those values must be provided.

* `rna.kallisto_fragment_length` Is the average fragment length.
* `rna.kallisto_sd_of_fragment_length` Is the standard deviation of the fragment lengths.

## Outputs

`Cromwell`: `Cromwell` will store outputs for each task under directory `cromwell-executions/[WORKFLOW_ID]/call-[TASK_NAME]/shard-[IDX]`. For all tasks `[IDX]` means a zero-based index for each replicate. In addition to the actual pipeline outputs, these directories contain a plethora of operational Cromwell-specific files. Use [croo](https://github.com/ENCODE-DCC/croo) to find and organize the pipeline outputs.

### Output files

#### Task Align

* `genomebam`, file name matches `*_genome.bam`. Bam aligned to genome.
* `annobam`, file name matches `*_anno.bam`. Bam aligned to annotation.
* `genome_flagstat` file name matches `*_genome_flagstat.txt`. Samtools flagstats on the genome bam.
* `genome_flagstat_json` file name matches `*_genome_flagstat.json`. Samtools flagstats on the genome bam in json format.
* `anno_flagstat` file name matches `*_anno_flagstat.txt`. Samtools flagstats on anno bam.
* `anno_flagstat_json` file name matches `*_anno_flagstat.json`. Samtools flagstats on anno bam in json format.
* `log` file name matches `*_Log.final.out`. STAR run log.
* `log_json` file name matches `*_Log.final.out`. STAR run log in json format.
* `python_log` file name is `align.log`. This file contains possible additional information on the pipeline step.

#### Task Kallisto

* `quants`, file name matches `*_abundance.tsv`. Kallisto quantifications.
* `python_log` file name is `kallisto_quant.log`. This file contains possible additional information on the pipeline step.

#### Task Bam to Signals

In case of an stranded run, the plus and minus strand signal tracks are separated (there will be four tracks per replicate).

* `unique`, file name matches `*niq.bw`. Contains the signal track of the uniquely mapped reads.
* `all`, the file name matches `*ll.bw`. Contains the signal track of all reads.
* `python_log` file name is `bam_to_signals.log`. This file contains possible additional information on the pipeline step.

#### Task RSEM Quant

* `genes_results`, file name matches `*.genes.results`. Contains gene quantifications.
* `isoforms_results`, file name matches `*.isoforms.results`. Contains isoform quantifications.
* `number_of_genes`, file name matches `*_number_of_genes_detected.json`. Contains the number of genes detected, which is determined as `TPM` value being greater than `1`.
* `python_log` file name is `rsem_quant.log`. This file contains possible additional information on the pipeline step.

#### Task Mad QC

This step is run if and only if the number of replicates is 2.

* `madQCplot`, file name matches `*_mad_plot.png`. Contains the MAD QC plot.
* `madQCmetrics` file name matches `*_mad_qc_metrics.json`. Contains MAD QC metrics.
* `python_log` file name is `mad_qc.log`. This file contains possible additional information on the pipeline step.

#### Task RNA QC

This step calculates additional metrics. At this time the only metric is to calculate reads by gene type. It is very **IMPORTANT** to look at the `python_log` of this step to see that the transcriptome bam did not contain any transcripts that are not present in the transcript ID to gene type mapping tsv. In case that happens, make sure you are using the STAR aligner and RSEM quantifier indexes you think you are using, and that all the other references are correct!

* `rnaQC`, file name matches `*_qc.json`. Contains additional QC metrics. For now the reads by gene type.
* `python_log` file name is `rna_qc.log`. This file contains **IMPORTANT** information on the pipeline step.

#### Note about resources:

The hardware resources needed to run the pipeline depend on the sequencing depth so it is hard to give definitive values that will be good for everyone. However, for every pipeline run, alignment is going to be the most memory-intensive task, quantitation with RSEM is going to be computationally hardest, kallisto will require some non-trivial amount of resources, and typically the rest of the tasks are rather light both in CPU and memory use. Disk usage again depends on the sequencing depth, but `"local-disk 100 HDD"` is a good starting point for all the tasks. The following are recommendations that are a sensible starting point for further optimizations in a typical case (non-CPU or memory related inputs omitted):

```
{
    "rna.align_ncpus" : 16,
    "rna.align_ramGB" : 60,
    "rna.rsem_ncpus" : 16,
    "rna.rsem_ramGB" : 60,
    "rna.kallisto_number_of_threads" : 8,
    "rna.kallisto_ramGB" : 30,
    "rna.bam_to_signals_ncpus" : 8,
    "rna.bam_to_signals_ramGB" : 30
}
```

In case of building index files for STAR and RSEM the sufficient amount of memory for GRCh38 is 60GB. The merge annotation workflow with gencode V29 annotation and tRNAs works with the default resources defined in the `merge_anno.wdl`. Smaller references may be able to run with less, but because index building typically needs to be done only once, it is most prudent to overshoot rather than waste time on several attempts. Following inputs with non-subsampled, full sized inputs have been tested on Google Cloud are good landmarks for defining the resources on other platforms as well:

Merge Annotation:

```
{
    "merge_anno.annotation" : "https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/@@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz",
    "merge_anno.tRNA" : "https://www.encodeproject.org/files/gencode.v29.tRNAs/@@download/gencode.v29.tRNAs.gtf.gz",
    "merge_anno.spikeins" : "https://github.com/ENCODE-DCC/rna-seq-pipeline/raw/master/test_data/ERCC_phiX.fa.gz",
    "merge_anno.output_filename" : "merged_annotation_V29.gtf.gz"
}
```

STAR index:

```
{
  "build_index.reference_sequence" : "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
  "build_index.spikeins" : "ERCC_phiX.fa.gz",
  "build_index.annotation" : "merged_annotation_V29.gtf.gz",
  "build_index.anno_version" : "v29",
  "build_index.genome" : "GRCh38",
  "build_index.index_type" : "prep_star",
  "build_index.ncpu" : 16,
  "build_index.memGB" : 60
}
```

RSEM index:

```
{
  "build_index.reference_sequence" : "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
  "build_index.spikeins" : "ERCC_phiX.fa.gz",
  "build_index.annotation" : "merged_annotation_V29.gtf.gz",
  "build_index.anno_version" : "v29",
  "build_index.genome" : "GRCh38",
  "build_index.index_type" : "prep_rsem",
  "build_index.ncpu" : 16,
  "build_index.memGB" : 60
}
```

Kallisto index:
```
{
    "build_index.reference_sequence" : "gencode.v29.transcripts_ERCC_phiX.fa.gz",
    "build_index.index_type" : "prep_kallisto",
    "build_index.ncpu" : 16,
    "build_index.memGB" : 60
}
```

```
{
  "rna.endedness": "single",
  "rna.fastqs_R1": [
    [
      "https://www.encodeproject.org/files/ENCFF824LLV/@@download/ENCFF824LLV.fastq.gz",
      "https://www.encodeproject.org/files/ENCFF729YAX/@@download/ENCFF729YAX.fastq.gz"
    ],
    [
      "https://www.encodeproject.org/files/ENCFF481BWJ/@@download/ENCFF481BWJ.fastq.gz",
      "https://www.encodeproject.org/files/ENCFF974EKR/@@download/ENCFF974EKR.fastq.gz"
    ]
  ],
  "rna.aligner": "star",
  "rna.bamroot": "ENCSR843RJV_testing",
  "rna.align_index": "gs://long_read_rna_production/caper_out/build_index/9320f158-6b7f-400a-899e-5b3f7e2ca46e/call-make_index/glob-08af9cb4bb3f985849894cfeafd3be7b/GRCh38_v29_ERCC_phiX_starIndex.tgz",
  "rna.rsem_index": "gs://long_read_rna_production/caper_out/build_index/b15d94e9-b4ef-49dc-b164-45000a0b3cb0/call-make_index/glob-08af9cb4bb3f985849894cfeafd3be7b/GRCh38_v29_ERCC_phiX_rsemIndex.tgz",
  "rna.strandedness": "unstranded",
  "rna.strandedness_direction": "unstranded",
  "rna.chrom_sizes": "gs://release-tests/Inputs/reference/GRCh38_EBV.chrom.sizes",
  "rna.align_ncpus": 16,
  "rna.align_ramGB": 60,
  "rna.kallisto_index": "gs://long_read_rna_production/caper_out/build_index/35ee9e2c-370d-416b-b08d-05307a63b7d1/call-make_index/glob-0a1c08388fcf0d2e32fc9d921d5a8910/gencode.v29.transcripts_ERCC_phiX.idx",
  "rna.kallisto_number_of_threads": 8,
  "rna.kallisto_ramGB": 30,
  "rna.kallisto_fragment_length": 250,
  "rna.kallisto_sd_of_fragment_length": 10,
  "rna.rna_qc_tr_id_to_gene_type_tsv": "gs://release-tests/Inputs/reference/gencodeV29pri-tRNAs-ERCC-phiX.transcript_id_to_genes.tsv",
  "rna.bam_to_signals_ncpus": 8,
  "rna.bam_to_signals_ramGB": 30,
  "rna.rsem_ncpus": 16,
  "rna.rsem_ramGB": 60,
  "rna.align_disk": "local-disk 200 SSD",
  "rna.kallisto_disk": "local-disk 200 SSD",
  "rna.rna_qc_disk": "local-disk 200 SSD",
  "rna.bam_to_signals_disk": "local-disk 200 SSD",
  "rna.mad_qc_disk": "local-disk 200 SSD",
  "rna.rsem_disk": "local-disk 200 SSD"
}
```
