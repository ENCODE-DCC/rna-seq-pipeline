#Dockerfile for ENCODE-DCC rna-seq-pipeline
#inspired by Jin Lee's atac-seq-pipeline docker image structure
FROM ubuntu:16.04
MAINTAINER Otto Jolanki 

RUN apt-get update
RUN apt-get install -y \
    python3-dev \
    python3-pip \
    wget

#Stick to Jin's way of organizing the directory structure
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

#additional installing here

RUN mkdir -p rna-seq-pipeline/src
COPY /src rna-seq-pipeline/src
ENV PATH="/software/rna-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]