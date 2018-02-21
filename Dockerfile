#Dockerfile for ENCODE-DCC rna-seq-pipeline
FROM ubuntu:16.04
MAINTAINER Otto Jolanki 

RUN apt-get update
RUN apt-get install -y \
    python3-dev \
    python3-pip \
    wget \
    git \
#libcurses is a samtools dependency
    libncurses5-dev \ 
#libboost is a tophat dependency
    libboost-all-dev


#Stick to Jin's way of organizing the directory structure
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

#Install STAR/Samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install
RUN wget http://bzip.org/1.0.6/bzip2-1.0.6.tar.gz && tar -xvf bzip2-1.0.6.tar.gz
RUN cd bzip2-1.0.6 && make && make install
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install

#Install STAR 2.5.1b
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

#Install Samtools 0.1.19
RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.9/samtools-0.1.9.tar.bz2 && tar -xvjf samtools-0.1.9.tar.bz2
RUN cd samtools-0.1.9 && make 
ENV PATH="/software/samtools-0.1.9:${PATH}"

RUN mkdir -p rna-seq-pipeline/src
#add a mount target dir for interactive testing
RUN mkdir mount_target
COPY /src rna-seq-pipeline/src
ENV PATH="/software/rna-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]