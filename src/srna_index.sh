#!/bin/bash -e

if [ $# -lt 4 ] || [ $# -gt 5 ]; then
    echo "usage v1: srna_index.sh <ref_fasta_gz> <annotation_gtf_gz> <annotation_version> <genome> [<gender>]"
    echo "Creates STAR index of reference genome and annotation for small-RNA-seq. Is independent of DX and encodeD."
    exit -1; 
fi
ref_fasta_gz=$1 # Reference genome assembly in gzipped fasta format.
anno_gtf_gz=$2  # Gene annotation in gzipped gtf format
anno=$3         # Annotation (e.g. 'v24')
genome=$4       # Genome (e.g. 'GRCh38')
if [ $# -eq 5 ]; then
    gender=$5        # Gender. Values: 'female', 'male', 'XX', 'XY' will be included in names.  Otherwise, gender neutral.  
fi

archive_file="${genome}_${anno}_sRNA_starIndex.tgz"
if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
    archive_file="${genome}_${gender}_${anno}_sRNA_starIndex.tgz"
fi
echo "-- Results will be: '${archive_file}'."

echo "-- Unzipping reference files..."
ref_fasta=${ref_fasta_gz%.gz}
ref_root=${ref_fasta%.fasta}
ref_root=${ref_root%.fa}
gunzip $ref_fasta_gz
anno_gtf=${anno_gtf_gz%.gz}
gunzip $anno_gtf_gz

echo "-- Build index for '${genome} and annotation ${anno}'..."
set -x
mkdir out
STAR --runMode genomeGenerate --genomeFastaFiles $ref_fasta --sjdbGTFfile $anno_gtf \
        --sjdbOverhang 1 --runThreadN 8 --genomeDir out/ --outFileNamePrefix out
set +x
    
# Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
echo "-- Create bam header..."
set -x
refComment="@CO\tREFID:$(basename ${ref_root})"
echo -e ${refComment} > out/star_bamCommentLines.txt
echo `cat "out/star_bamCommentLines.txt"`
set +x

echo "-- Create archive file..."
set -x
tar -czvf $archive_file out/
set +x

echo "-- The results..."
ls -l $archive_file

