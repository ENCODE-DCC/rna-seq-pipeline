# INSTALLATION

To run the pipeline you need to install following software. Running the pipeline on Google Cloud or DNA Nexus requires additional setup detailed below.

## Java 8

Java is required to run execution engine [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution) and WDL-to-DNANexus compiler [dxWDL](https://github.com/dnanexus/dxWDL).
To check which Java version you already have, run:
```bash
  $ java -version
```
You are looking for 1.8 or higher. If the requirement is not fulfilled follow installation instructions for [mac](https://java.com/en/download/help/mac_install.xml) or
[linux](http://openjdk.java.net/install/) or use your favorite installation method.

## Docker

Pipeline code is packaged and distributed in Docker containers, and thus Docker installation is needed. 
Follow instructions for [mac](https://docs.docker.com/docker-for-mac/install/) or [linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/#upgrade-docker-after-using-the-convenience-script).

##Singularity

If for some reason (most likely because of limitations on your HPC resource) you cannot run Docker, install [singularity](https://www.sylabs.io/guides/2.6/user-guide/installation.html) and have a look at [HOWTO](howto.md#local-with-singularity) for an example of how to run pipeline with singularity. Pipeline requires singularity version `>=2.5.2`, the link takes you to version `2.6.0`.

## Google Cloud

If you are intending to run the pipeline on Google Cloud platform, the following setup is needed:

1. Sign up for a Google account.
2. Go to [Google Project](https://console.developers.google.com/project) page and click "SIGN UP FOR FREE TRIAL" on the top left and agree to terms.
3. Set up a payment method and click "START MY FREE TRIAL".
4. Create a [Google Project](https://console.developers.google.com/project) `[YOUR_PROJECT_NAME]` and choose it on the top of the page.
5. Create a [Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser) `gs://[YOUR_BUCKET_NAME]` by clicking on a button "CREATE BUCKET" and create it to store pipeline outputs.
6. Find and enable following APIs in your [API Manager](https://console.developers.google.com/apis/library). Click a back button on your web brower after enabling each.
    * Compute Engine API
    * Google Cloud Storage (DO NOT click on "Create credentials")
    * Google Cloud Storage JSON API
    * Genomics API

7. Install [Google Cloud Platform SDK](https://cloud.google.com/sdk/downloads) and authenticate through it. You will be asked to enter verification keys. Get keys from the URLs they provide.
    ```
      $ gcloud auth login --no-launch-browser
      $ gcloud auth application-default login --no-launch-browser
    ```

8. If you see permission errors at runtime, then unset environment variable `GOOGLE_APPLICATION_CREDENTIALS` or add it to your BASH startup scripts (`$HOME/.bashrc` or `$HOME/.bash_profile`).
    ```
      unset GOOGLE_APPLICATION_CREDENTIALS
    ```

9. Set your default Google Cloud Project. Pipeline will provision instances on this project.
    ```
      $ gcloud config set project [YOUR_PROJECT_NAME]
    ```

* For an example on how to run the pipeline on Google Cloud, see [HOWTO](howto.md#google-cloud).

## DNA Nexus

If you are intending to build and run the pipeline on DNA Nexus, the following setup is needed:

1. Sign up for a [DNANexus account](https://platform.dnanexus.com/register).

2. Create a new [DX project](https://platform.dnanexus.com/projects) with name `[YOUR_PROJECT_NAME]` by clicking on "+New Project" on the top left.

3. Install [DNA Nexus SDK](https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK). Run `dx login` to log in to the platform using your credentials.

4. Download dxWDL:
    ```bash 
      $ wget https://github.com/dnanexus/dxWDL/releases/download/0.75/dxWDL-0.75.jar
    ```

5. Make dxWDL executable:
    ```bash
      $ chmod 755 dxWDL-0.75.jar
    ```

* For an example on how to run the pipeline on DNA Nexus, see [HOWTO](howto.md#dna-nexus)