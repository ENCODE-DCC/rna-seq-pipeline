# INSTALLATION

To run the pipeline you need to do some setup. The exact steps you need to take depend on the platform you are running the pipeline on, and will be detailed below and in [HOWTO](howto.md). Independent of platform, running the pipeline is done using [caper](https://github.com/ENCODE-DCC/caper) and (optional but recommended) output organization is done using [croo](https://github.com/ENCODE-DCC/croo). Both `caper` and `croo` require `python` version 3.4.1 or newer.

## Caper

Direct usage of the execution engine [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution) features complicated backend configuration, workflow options and command line parameters. Caper hides the complexity and consolidates configuration in one file. Caper is available in [PyPI](https://pypi.org/project/caper/) and it is installed by running:

```bash
  $ pip install caper
```

Note that conda run mode that is described in caper documentation is not supported by this pipeline.

## Croo

The way [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution) organizes pipeline outputs is not always the most clear and findable. Croo is a tool to reorganize the files in more readable manner. Croo is available in [PyPI](https://pypi.org/project/croo/) and it is installed by running:

```bash
  $ pip install croo
```

## Java 8

Java is required to run execution engine [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution) that `caper` uses under the hood.
To check which Java version you already have, run:
```bash
  $ java -version
```
You are looking for 1.8 or higher. If the requirement is not fulfilled follow installation instructions for [mac](https://java.com/en/download/help/mac_install.xml) or
[linux](http://openjdk.java.net/install/) or use your favorite installation method.

## Docker

Pipeline code is packaged and distributed in Docker containers, and thus Docker installation is needed.
Follow instructions for [mac](https://docs.docker.com/docker-for-mac/install/) or [linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/#upgrade-docker-after-using-the-convenience-script).

## Singularity

If you want to use Singularity instead of Docker, install [singularity](https://www.sylabs.io/guides/3.1/user-guide/installation.html). Pipeline requires singularity version `>=2.5.2`, the link takes you to version `3.1`.

## Google Cloud

If you are intending to run the pipeline on Google Cloud platform, follow the [caper setup instructions for GCP](https://github.com/ENCODE-DCC/caper/blob/master/docs/conf_gcp.md).
* For an example on how to run the pipeline on Google Cloud, see [HOWTO](howto.md#google-cloud).

## AWS

If you are intending to run the pipeline on AWS, follow the [caper setup instructions for AWS](https://github.com/ENCODE-DCC/caper/blob/master/docs/conf_aws.md).

## Cromwell (optional)

We recommend using `caper` for running the pipeline, although it is possible to use Cromwell directly. Backend file and workflow options files necessary for direct Cromwell use are included in the repository for local testing purposes, but they are not actively maintained to follow cloud API changes etc.
