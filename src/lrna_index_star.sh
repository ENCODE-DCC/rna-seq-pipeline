#!/bin/bash -e

if [ $# -lt 5 ] || [ $# -gt 6 ]; then
    echo "usage v1: lrna_index_star.sh <ref_fasta_gz> <spike_in_fasta_gz> <annotation_gtf_gz> <annotation_version> <genome> [<gender>]"
    echo "Creates STAR index of reference genome, annotation and spike-ins for long-RNA-seq. Is independent of DX and encodeD."
    exit -1; 
fi
ref_fasta_gz=$1      # Reference genome assembly in gzipped fasta format.
spike_fasta_gz=$2    # All spike-ins in single gzipped fasta format.
anno_gtf_gz=$3       # Gene annotation in gzipped gtf format
anno=$4              # Annotation (e.g. 'v24')
genome=$5            # Genome (e.g. 'GRCh38')
if [ $# -eq 6 ]; then
    gender=$6        # Gender. Values: 'female', 'male', 'XX', 'XY' will be included in names.  Otherwise, gender neutral.  
fi

echo "-- Unzipping reference files..."
ref_fasta=${ref_fasta_gz%.gz}
ref_root=${ref_fasta%.fasta}
ref_root=${ref_root%.fa}
gunzip $ref_fasta_gz
spike_fasta=${spike_fasta_gz%.gz}
spike_root=${spike_fasta%.fasta}
spike_root=${spike_root%.fa}
gunzip $spike_fasta_gz
anno_gtf=${anno_gtf_gz%.gz}
anno_root=${anno_gtf%.gtf}
gunzip $anno_gtf_gz

archive_file="${genome}_${anno}_${spike_root}_starIndex.tgz"
if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
    archive_file="${genome}_${gender}_${anno}_${spike_root}_starIndex.tgz"
fi
echo "-- Results will be: '${archive_file}'."

echo "-- Build index..."
set -x
mkdir out
STAR --runMode genomeGenerate --genomeFastaFiles $ref_fasta $spike_fasta \
     --sjdbOverhang 100 --sjdbGTFfile $anno_gtf --runThreadN 8 --genomeDir out/ \
     --outFileNamePrefix out
set +x

# Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
echo "-- Create bam header..."
set -x
refComment="@CO\tREFID:$(basename ${ref_root})"
annotationComment="@CO\tANNID:$(basename ${anno_root})"
spikeInComment="@CO\tSPIKEID:${spike_root}"
echo -e ${refComment} > out/star_bamCommentLines.txt
echo -e ${annotationComment} >> out/star_bamCommentLines.txt
echo -e ${spikeInComment} >> out/star_bamCommentLines.txt
echo `cat "out/star_bamCommentLines.txt"`
set +x

echo "-- Tar up results..."
set -x
tar -czvf $archive_file out/
set +x

echo "-- The results..."
ls -l $archive_file

