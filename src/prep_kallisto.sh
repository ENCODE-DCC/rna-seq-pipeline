#!/bin/bash -e

ref_fasta_gz=$1
ref_base = $(basename ref_fasta_gz)
index_root=${ref_base%.fa*.gz}.idx

kallisto index -i $index_root $ref_fasta_gz
