# HOWTO

Here are recipes for running analyses on different platforms.
Before following these instructions, make sure you have completed installation and possible account setup detailed in [installation instructions](installation.md). Note that although running the pipeline directly with Cromwell is still possible, using [caper](https://github.com/ENCODE-DCC/caper) is the canonical, supported and official way to use ENCODE Uniform Processing Pipelines. The examples below use command `caper run`, which is the simplest way to run a single pipeline instance. For running multiple pipelines in production setting we recommend using caper server. To find details on setting up the server, refer to [caper documentation](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#usage).

Note that the files used in these exampled are first restricted to reads from chromosome 19, and then further subsampled to 10000 reads. The cpu and memory resources reflect the size of inputs. For resource guidelines with full sized data, see discussion [here](reference.md#note-about-resources).


# CONTENTS

[Google Cloud](howto.md#google-cloud)  
[Local with Docker](howto.md#local-with-docker)  
[Local with Singularity](howto.md#local-with-singularity)  
[Other Platforms](howto.md#other-platforms)  
[Building index files](howto.md#building-indexes)  

# RUNNING THE PIPELINE

## Google Cloud

The goal is to run a paired-end, strand-specific experiment on Google Cloud Platform.
Make sure you have completed the steps for installation and Google Cloud setup described in the [installation instructions](installation.md#google-cloud). The following assumes your Google Cloud project is `[YOUR_PROJECT]`, you have created a bucket `gs://[YOUR_BUCKET_NAME]`, and also directories `inputs`, `output` and `reference` in the bucket.

1. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

2. Get STAR and kallisto index files:

```bash
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz -o test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx -o test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx
```

3. Copy indexes and input data into the cloud:

```bash
  $ gsutil cp test_data/ENCSR653DFZ* gs://[YOUR_BUCKET_NAME]/inputs/
  $ gsutil cp test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz gs://[YOUR_BUCKET_NAME]/reference/
  $ gsutil cp test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx gs://[YOUR_BUCKET_NAME]/reference/
  $ gsutil cp test_data/GRCh38_v24_ERCC_phiX_rsemIndex_chr19only.tgz gs://[YOUR_BUCKET_NAME]/reference/
  $ gsutil cp test_data/GRCh38_EBV.chrom.sizes gs://[YOUR_BUCKET_NAME]/reference/
```

4. Set up the `input.json`:

    Copy the following into `input.json` in your favorite text editor.

```
{
    "rna.endedness" : "paired",
    "rna.fastqs_R1" : ["gs://[YOUR_BUCKET_NAME]/inputs/ENCSR653DFZ_rep1_chr19_10000reads_R1.fastq.gz", "gs://[YOUR_BUCKET_NAME]/inputs/ENCSR653DFZ_rep2_chr19_10000reads_R1.fastq.gz"],
    "rna.fastqs_R2" : ["gs://[YOUR_BUCKET_NAME]/inputs/ENCSR653DFZ_rep1_chr19_10000reads_R2.fastq.gz", "gs://[YOUR_BUCKET_NAME]/inputs/ENCSR653DFZ_rep2_chr19_10000reads_R2.fastq.gz"],
    "rna.aligner" : "star",
    "rna.align_index" : "gs://[YOUR_BUCKET_NAME]/reference/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz",
    "rna.rsem_index" : "gs://[YOUR_BUCKET_NAME]/reference/GRCh38_v24_ERCC_phiX_rsemIndex_chr19only.tgz",
    "rna.kallisto_index" : "gs://[YOUR_BUCKET_NAME]/reference/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx",
    "rna.bamroot" : "PE_stranded",
    "rna.strandedness" : "stranded",
    "rna.strandedness_direction" : "reverse",
    "rna.chrom_sizes" : "gs://[YOUR_BUCKET_NAME]/reference/GRCh38_EBV.chrom.sizes",
    "rna.align_ncpus" : 2,
    "rna.align_ramGB" : 4,
    "rna.rsem_ncpus" : 2,
    "rna.rsem_ramGB" : 4,
    "rna.kallisto_number_of_threads" : 2,
    "rna.kallisto_ramGB" : 4,
    "rna.rna_qc_tr_id_to_gene_type_tsv" : "gs://[YOUR_BUCKET_NAME]/reference/gencodeV24pri-tRNAs-ERCC-phiX.transcript_id_to_genes.tsv",
    "rna.bam_to_signals_ncpus" : 1,
    "rna.bam_to_signals_ramGB" : 2,
    "rna.align_disk" : "local-disk 20 HDD",
    "rna.kallisto_disk" : "local-disk 20 HDD",
    "rna.rna_qc_disk" : "local-disk 20 HDD",
    "rna.bam_to_signals_disk" : "local-disk 20 HDD",
    "rna.mad_qc_disk" : "local-disk 20 HDD",
    "rna.rsem_disk" : "local-disk 20 HDD"
}
```

Replace `[YOUR_BUCKET_NAME]` with the name of the bucket you created.

5. Run the pipeline using caper:

```bash
  $ caper run rna-seq-pipeline.wdl -i input.json -b gcp -m testrun_metadata.json
```

6. Run croo, to to make finding outputs easier:

```bash
  $ croo testrun_metadata.json --out-dir gs://[YOUR_BUCKET_NAME]/croo_out
```

This command will output into the bucket an HTML table, that shows the locations of the outputs nicely organized. Note that if your output bucket is not public, you need to be logged into your google account to be able to follow the links.


## Local with Docker

The goal is to run a single-end, non-strand-specific experiment on a local computer.

1. Get the code:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

2. Get STAR and kallisto index files:

```bash
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz -o test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx -o test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx
```

The other data that is required to complete this recipe is included in the repository within test_data directory.

3. Set up the `input.json`:

    Copy the following into `input.json` in your favorite text editor.

```
{
    "rna.endedness" : "single",
    "rna.fastqs_R1" : ["[PATH_TO_REPO]/rna-seq-pipeline/test_data/rep1_ENCSR510QZW_chr19only_10000_reads.fastq.gz","<path-to-repo>/rna-seq-pipeline/test_data/rep2_ENCSR510QZW_chr19only_10000_reads.fastq.gz"],
    "rna.aligner" : "star",
    "rna.bamroot" : "SE_unstranded",
    "rna.align_index" : "[PATH_TO_REPO]/rna-seq-pipeline/test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz",
    "rna.rsem_index" : "[PATH_TO_REPO]/rna-seq-pipeline/test_data/GRCh38_v24_ERCC_phiX_rsemIndex_chr19only.tgz",
    "rna.kallisto_index" : "[PATH_TO_REPO]/rna-seq-pipeline/test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx",
    "rna.strandedness" : "unstranded",
    "rna.strandedness_direction" : "unstranded",
    "rna.chrom_sizes" : "[PATH_TO_REPO]/rna-seq-pipeline/test_data/GRCh38_EBV.chrom.sizes",
    "rna.align_ncpus" : 2,
    "rna.align_ramGB" : 4,
    "rna.rsem_ncpus" : 2,
    "rna.rsem_ramGB" : 4,
    "rna.kallisto_number_of_threads" : 2,
    "rna.kallisto_ramGB" : 4,
    "rna.kallisto_fragment_length" : 250,
    "rna.kallisto_sd_of_fragment_length" : 10,
    "rna.rna_qc_tr_id_to_gene_type_tsv" : "[PATH_TO_REPO]/rna-seq-pipeline/transcript_id_to_gene_type_mappings/gencodeV24pri-tRNAs-ERCC-phiX.transcript_id_to_genes.tsv",
    "rna.bam_to_signals_ncpus" : 1,
    "rna.bam_to_signals_ramGB" : 2,
    "rna.align_disk" : "local-disk 20 HDD",
    "rna.kallisto_disk" : "local-disk 20 HDD",
    "rna.rna_qc_disk" : "local-disk 20 HDD",
    "rna.bam_to_signals_disk" : "local-disk 20 HDD",
    "rna.mad_qc_disk" : "local-disk 20 HDD",
    "rna.rsem_disk" : "local-disk 20 HDD"

}
```

Replace `[PATH_TO_REPO]` with the location you cloned the code into.


4. Run the pipeline using caper:

```
  $ caper run rna-seq-pipeline.wdl -i input.json -m testrun_metadata.json --docker
```

5. Organize outputs with croo:

```bash
  $ croo testrun_metadata.json --out-dir organized_outputs
```

This will create directory tree grouped by task and replicate. Disk space does not get wasted, by default croo creates symbolic links in local mode.

## Local with Singularity

The goal is to run a single-end non-strand-specific experiment locally using singularity.

1. Make sure you have singularity version greater or equal to `2.5.2` installed in your system.

2. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

3. Get STAR and kallisto index files:

```bash
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz -o test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz
  $ curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx -o test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix_k31_kallisto.idx
```

4. Run the pipeline using caper:
```
$ caper run rna-seq-pipeline.wdl -i test/test_workflow/SE_unstranded_input.json --singularity
```

5. Organize outputs with croo:

```bash
  $ croo testrun_metadata.json --out-dir organized_outputs
```

This will create directory tree grouped by task and replicate. Disk space does not get wasted, by default croo creates symbolic links in local mode.


# Other platforms

Running on other platforms is similar, because the caper takes care of the details for you. See [caper documentation](https://github.com/ENCODE-DCC/caper#installation) for further details.

# BUILDING INDEXES

Most likely you will not need to build STAR, RSEM or Kallisto indexes if you are working with data from human or mouse samples and are using standard spikeins. Links to available reference files can be found [here](reference.md#genome-reference-files). If you need to build references from scratch, you can find the needed workflow code in `per_task_wdl` directory in this repo. See [example inputs](reference.md#note-about-resources) for guidance.
