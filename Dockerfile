# Dockerfile for ENCODE-DCC rna-seq-pipeline
FROM ubuntu:16.04
MAINTAINER Otto Jolanki 

RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-pip \
    wget \
    git \
    unzip \
    bzip2 \
# libcurses is a samtools dependency
    libncurses5-dev \ 
    r-base-core \
    ghostscript && rm -rf /var/lib/apt/lists/*

# Add mountpoint directory
RUN mkdir /data
# Stick to Jin's way of organizing the directory structure
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install STAR/Samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install && rm ../zlib-1.2.11.tar.gz

RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install && rm ../xz-5.2.3.tar.gz

# Install STAR 2.5.1b
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR && rm ../2.5.1b.tar.gz
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

# Install Kallisto 0.44.0
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && tar -xzf kallisto_linux-v0.44.0.tar.gz
ENV PATH="/software/kallisto_linux-v0.44.0:${PATH}"

# Install Samtools 0.1.19
RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 && tar -xvjf samtools-0.1.19.tar.bz2
# Bzip unloads the directory with privs 750 instead of 755. This trips singularity that runs as an unprivileged user.
RUN chmod 755 samtools-0.1.19
RUN cd samtools-0.1.19 && make && rm ../samtools-0.1.19.tar.bz2
ENV PATH="/software/samtools-0.1.19:${PATH}"

# Install RSEM 1.2.31
RUN wget https://github.com/deweylab/RSEM/archive/v1.2.31.zip
RUN unzip v1.2.31.zip && rm v1.2.31.zip
RUN cd RSEM-1.2.31 && make
ENV PATH="/software/RSEM-1.2.31:${PATH}"

# Install BedGraphToBigWig and bedSort
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && chmod +x bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort && chmod +x bedSort

RUN mkdir -p rna-seq-pipeline/src
COPY /src rna-seq-pipeline/src
COPY logger_config.json /software/
ENV PATH="/software/rna-seq-pipeline/src:${PATH}"
ARG GIT_COMMIT_HASH
ENV GIT_HASH=${GIT_COMMIT_HASH}
ARG BRANCH
ENV BUILD_BRANCH=${BRANCH}
ARG BUILD_TAG
ENV MY_TAG=${BUILD_TAG}

ENTRYPOINT ["/bin/bash", "-c"]