# HOWTO

Here are recipes for running analyses on different platforms.
Before following these instructions, make sure you have completed installation and possible account setup detailed in [installation instructions](installation.md). Note that although running the pipeline directly with Cromwell is still possible, using [caper](https://github.com/ENCODE-DCC/caper) is the canonical, supported and official way to use ENCODE Uniform Processing Pipelines. The examples below use command `caper run`, which is the simplest way to run a single pipeline instance. For running multiple pipelines in production setting we recommend using caper server. To find details on setting up the server, refer to [caper documentation](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#usage).

Note that the files used in these exampled are first restricted to reads from chromosome 19, and then further subsampled to 10000 reads. The cpu and memory resources reflect the size of inputs. For resource guidelines with full sized data, see discussion [here](reference.md#note-about-resources).


# CONTENTS

## Running Analyses

[Google Cloud](howto.md#google-cloud)  
[Local with Docker](howto.md#local-with-docker)  
[Local with Singularity](howto.md#local-with-singularity)  
[Sherlock with Singularity](howto.md#sherlock-with-singularity)  
[SLURM](howto.md#slurm)  

## Building indexes

[Merge Annotation](howto.md#merge-annotation)  
[Build STAR Index](howto.md#build-star-index)  
[Build RSEM Index](howto.md#build-rsem-index)  
[Build Kallisto Index](howto.md#build-kallisto-index)  

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

6. Run croo, to easier view the pipeline output:

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

5. See outputs in `cromwell-executions/rna/[RUNHASH]`.


## Sherlock with Singularity

The goal is to run a paired-end, strand-specific experiment on Sherlock using singularity.

1. SSH into Sherlock's login node:

```bash
  $ ssh user@login.sherlock.stanford.edu
```

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

4. Load singularity and java modules into your environment. You can add the lines to your `~/.bashrc` or `~/.bash_profile` if you want them available always when you log in to Sherlock:

```bash
  $ module load system singularity
  $ module load java
```

5. Build the singularity image for the pipeline. The following pulls the pipeline docker image, and uses that to construct the singularity image. The image will be stored in `~/.singularity`. It is bad practice to build images (or do any other intensive work) on login nodes. For this reason we will first invoke an interactive session on a different node by running `sdev` command, and building there (It will take few seconds to get back into the shell after running `sdev`).

```bash
  $ sdev
  $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/rna-seq-pipeline:v1.0
  $ exit #this takes you back to the login node
```

6. Open `workflow_opts/sherlock.json` in your favorite text editor:

```
{
    "default_runtime_attributes" : {
        "singularity_container" : "~/.singularity/rna-seq-pipeline-v1.0.simg",
        "singularity_command_options" : "--bind /scratch,/lscratch,/oak/stanford"
    }
}
```

The default SLURM partition is `normal`. If you want to use some other partition, as you probably will when running a full sized experiment, add `"slurm_partition"` line in the workflow options. After this addition your `workflow_opts/sherlock.json` looks like this:

```
{
    "default_runtime_attributes" : {
        "singularity_container" : "~/.singularity/rna-seq-pipeline-v1.0.simg",
        "slurm_partition" : "SLURM_PARTITION"
        "singularity_command_options" : "--bind /scratch,/lscratch,/oak/stanford"
    }
}
```

7. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm_singularity cromwell-34.jar run rna-seq-pipeline.wdl -i test/test_workflow/PE_stranded_input.json -o workflow_opts/sherlock.json
```

8. See the outputs in `cromwell-executions/rna/[RUNHASH]`.

# SLURM

Using a generic SLURM cluster should be quite similar to the Stanford Sherlock as documented above (which is a SLURM machine of a specific kind). The main differences are that you may need to [install](installation.md#singularity) singularity and edit `workflow_opts/slurm.json` to include your information and the directories that contain your input data.

* Singularity version has to be `>=2.5.2`

The `slurm.json` template file looks like this:

```
{
    "default_runtime_attributes" : {
      "slurm_partition": "[YOUR_SLURM_PARTITION]",
      "slurm_account": "[YOUR_SLURM_ACCOUNT]",
      "singularity_container" : "~/.singularity/rna-seq-pipeline-v1.0.simg",
      "singularity_command_options" : "--bind /your/,[DATA_DIR1],[DATA_DIR2],..."
    }
}
```

1. Setup your partition, account and data:

Set your partition/account in workflow_opts/slurm.json. If your SLURM cluster does not require either user's partition or account information, then remove them from this file. Otherwise, `YOUR_SLURM_PARTITON` or `YOUR_SLURM_ACCOUNT` will be used internally for `srun ... --partition YOUR_SLURM_PARTITION` and `srun ... --account YOUR_SLURM_PARTITION`, respectively.

2. Build the singularity image:

```bash
  $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/rna-seq-pipeline:v1.0
```

3. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm_singularity cromwell-34.jar run rna-seq-pipeline.wdl -i [INPUT] -o workflow_opts/slurm.json
```


On the `"singularity_command_options"` line, add the paths to the directories that contain your input data. This is comma separated and can include several items.

# BUILDING INDEXES

Most likely you will not need to build STAR or RSEM indexes if you are working with data from human or mouse samples and are using standard spikeins. STAR indexes can be downloaded [here](https://www.encodeproject.org/references/ENCSR314WMD/) and RSEM indexes [here](https://www.encodeproject.org/references/ENCSR219BJA/). Genome references used in building the aforementioned indexes can be found [here](https://www.encodeproject.org/references/ENCSR425FOI/).

## Merge Annotation

This step precedes the STAR and RSEM index building steps. This step takes a gene annotation (e.g. gencode.v19.annotation.gtf.gz) gzipped gtf file, the corresponding tRNA gzipped gtf and a spike-in set (e.g. ERCC) gzipped fasta file and combines them into a single merged annotation gzipped gtf file. This file will be used as input to all three 'prep' indexing steps.

The goal is to run the Merge Annotation step on a local machine using Docker.

1. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

2. Find input file `merge_anno_input.json` in folder `input_json_templates/per_task_inputs`. The file looks like this:

```
{
    "merge_anno.annotation" : "test_data/gencode.v24.primary_assembly.annotation.gtf.gz",
    "merge_anno.tRNA" : "test_data/gencode.v24.tRNAs.gtf.gz",
    "merge_anno.spikeins" : "test_data/ERCC_phiX.fa.gz",
    "merge_anno.output_filename" : "merged_annotation.gtf.gz"
}
```
There is no need to edit this file.

3. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=Local cromwell-34.jar run per_task_wdl/merge_anno.wdl -i input_json_templates/per_task_inputs/merge_anno_input.json -o workflow_opts/docker.json
```

4. Find the output in `cromwell-executions/merge_anno/[RUNHASH]`

## Build the STAR Index

The goal is to build on the [previous step](howto.md#merge-annotation) and build a STAR index, restricted to chromosome 19, using a local machine with Docker.

1. Make sure you have run the [previous step](howto.md#merge-annotation) and have located the output (`merged_annotation.gtf.gz`) of that step.

2. Find the input file `build_index_STAR.json` in `input_json_templates/per_task_inputs`. The file looks like this:

```
{
  "build_index.reference_sequence" : "test_data/GRCh38_no_alt_analysis_set_GCA_000001405.15_onlychr19.fa.gz",
  "build_index.spikeins" : "test_data/ERCC_phiX.fa.gz",
  "build_index.annotation" : "[path-to-output]/gencodeV24pri-tRNAs-ERCC-phiX_onlychr19_and_spikeins.gtf.gz",
  "build_index.anno_version" : "v24",
  "build_index.genome" : "GRCh38",
  "build_index.index_type" : "prep_star"
}
```

Open the file in your favorite text editor, and replace `[path-to-output]` with the path to the merged annotation file produced in the [previous step](howto.md#merge-annotation).

3. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=Local cromwell-34.jar run per_task_wdl/build_genome_index.wdl -i input_json_templates/per_task_inputs/build_index_STAR.json -o workflow_opts/docker.json
```

4. Find the output in `cromwell-executions/build_genome_index/[RUNHASH]`

## Build the RSEM Index

The goal is to build on the [previous step](howto.md#merge-annotation) and build a RSEM index, restricted to chromosome 19, using a local machine with Docker.

1. Make sure you have run the [previous step](howto.md#merge-annotation) and have located the output (`merged_annotation.gtf.gz`) of that step.

2. Find the input file `build_index_RSEM.json` in `input_json_templates/per_task_inputs`. The file looks like this:

```
{
  "build_index.reference_sequence" : "test_data/GRCh38_no_alt_analysis_set_GCA_000001405.15_onlychr19.fa.gz",
  "build_index.spikeins" : "test_data/ERCC_phiX.fa.gz",
  "build_index.annotation" : "[path-to-output]/gencodeV24pri-tRNAs-ERCC-phiX_onlychr19_and_spikeins.gtf.gz",
  "build_index.anno_version" : "v24",
  "build_index.genome" : "GRCh38",
  "build_index.index_type" : "prep_rsem"
}
```

Open the file in your favorite text editor and replace `[path-to-output]` with the path to the merged annotation file produced in the [previous step](howto.md#merge-annotation).

3. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=Local cromwell-34.jar run per_task_wdl/build_genome_index.wdl -i input_json_templates/per_task_inputs/build_index_RSEM.json -o workflow_opts/docker.json
```

4. Find the output in `cromwell-executions/build_genome_index/[RUNHASH]`

## Build Kallisto Index

The goal is to build index for kallisto using a local machine with Docker.

1. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

2. Find the input file `build_index_Kallisto.json` in `input_json_templates/per_task_inputs`. The file looks like this:

```
{
    "build_index.reference_sequence" : "test_data/Homo_sapiens.GRCh38.cdna.all.chr19_ERCC_phix.fa.gz",
    "build_index.index_type" : "prep_kallisto",
}
```

As you can see, this file has fewer inputs than [STAR](howto.md#build-star-index) and [RSEM](howto.md#build-rsem-index) steps. The reason for this is that when building the index kallisto uses only the transcriptome and does not need annotations. Additionally, the spikein sequences are concatenated into the reference file, instead of providing them in a separate input. The inputs used in this example are [ERCC spikes](https://www.encodeproject.org/files/ENCFF001RTP/) and [Human cDNA](ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz), which is restricted to chromosome 19.

3. Run the pipeline:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=Local cromwell-34.jar run per_task_wdl/build_genome_index.wdl -i input_json_templates/per_task_inputs/build_index_Kallisto.json -o workflow_opts/docker.json
```

4. Find the output in `cromwell-executions/build_genome_index/[RUNHASH]`
