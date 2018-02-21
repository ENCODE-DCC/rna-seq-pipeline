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

#Install STAR dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install

#Install STAR 2.5.1b
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

RUN mkdir -p rna-seq-pipeline/src
#add a mount target dir for interactive testing
RUN mkdir mount_target
COPY /src rna-seq-pipeline/src
ENV PATH="/software/rna-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]