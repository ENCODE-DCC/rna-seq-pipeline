# HOWTO

Here are recipes for running analyses on different platforms.
Before following these instructions, make sure you have completed installation and possible account setup detailed in [installation instructions](installation.md). Note that although running the pipeline directly with Cromwell is still possible, using [caper](https://github.com/ENCODE-DCC/caper) is the canonical, supported and official way to use ENCODE Uniform Processing Pipelines. The examples below use command `caper run`, which is the simplest way to run a single pipeline instance. For running multiple pipelines in production setting we recommend using caper server. To find details on setting up the server, refer to [caper documentation](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#usage).

Note that the files used in these exampled are first restricted to reads from chromosome 19, and then further subsampled to 10000 reads. The cpu and memory resources reflect the size of inputs. For resource guidelines with full sized data, see discussion [here](reference.md#note-about-resources).


# CONTENTS

[Running the pipeline](howto.md#running-the-pipeline)
[Truwl](howto.md#truwl)
[Building index files](howto.md#building-indexes)  


# RUNNING THE PIPELINE

The goal is to run a paired-end, strand-specific experiment on a chosen platform. For GCP, make sure you have completed the steps for installation and Google Cloud setup described in the [installation instructions](installation.md#google-cloud). The following assumes your Google Cloud project is `[YOUR_PROJECT]`, you have created a bucket `gs://[YOUR_BUCKET_NAME]`

You can use this input JSON (`gs://encode-pipeline-test-samples/rna-seq-pipeline/input.json`) universally on all platforms. Caper will localize all files in the input JSON according to the chose platform.

1. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/rna-seq-pipeline
  $ cd rna-seq-pipeline
```

2. Run it on a chosen platform.

```bash
  # Run it on GCP
  $ caper run rna-seq-pipeline.wdl -m testrun_metadata.json -i gs://encode-pipeline-test-samples/rna-seq-pipeline/input.json -b gcp 

  # Run it locally with docker
  $ caper run rna-seq-pipeline.wdl -m testrun_metadata.json -i gs://encode-pipeline-test-samples/rna-seq-pipeline/input.json -b local --docker

  # Run it locally with singularity
  $ caper run rna-seq-pipeline.wdl -m testrun_metadata.json -i gs://encode-pipeline-test-samples/rna-seq-pipeline/input.json -b local --singularity

  # sbatch with singularity (read caper's doc very carefully)
  $ sbatch [...SBATCH PARTITON/ACCOUNT/RESOURCE PARAMS HERE...] --wrap "caper run rna-seq-pipeline.wdl -m testrun_metadata.json -i gs://encode-pipeline-test-samples/rna-seq-pipeline/input.json -b slurm --singularity"
```

3. Organize outputs with `croo`. This command will output into the bucket (or local directory) an HTML table, that shows the locations of the outputs nicely organized. For GCP, note that if your output bucket is not public, you need to be logged into your google account to be able to follow the links.

```bash
  # Organize outputs on GCP (this will not transfer/copy any original outputs)
  $ croo testrun_metadata gs://[YOUR_BUCKET_NAME]/croo_out

  # Organize outputs locally (this will not transfer/copy any original outputs, it will simply soft-link)
  $ croo testrun_metadata --out-dir organized_outputs
```

# Truwl

You can run this pipeline on [truwl.com](https://truwl.com/). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress in
a ready-to-go environment. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right per
missions. You can see the paired-end stranded example case from this repo [here](https://truwl.com/workflows/library/RNA-seq%20pipeline/v1.2.1/instances/WF_0741dc.ac.1e10). The example j
ob (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.j
son file.

# BUILDING INDEXES

Most likely you will not need to build STAR, RSEM or Kallisto indexes if you are working with data from human or mouse samples and are using standard spikeins. Links to available reference files can be found [here](reference.md#genome-reference-files). If you need to build references from scratch, you can find the needed workflow code in `per_task_wdl` directory in this repo. See [example inputs](reference.md#note-about-resources) for guidance.

