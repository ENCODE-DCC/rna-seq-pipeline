#!/bin/bash -e

if [ $# -lt 5 ] || [ $# -gt 6 ]; then
    echo "usage v1: lrna_index_rsem.sh <ref_fasta_gz> <spike_in_fasta_gz> <annotation_gtf_gz> <annotation_version> <genome> [<gender>]"
    echo "Creates RSEM index of reference genome, annotation and spike-ins for long-RNA-seq. Is independent of DX and encodeD."
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
gunzip -c $ref_fasta_gz > $ref_fasta
spike_fasta=${spike_fasta_gz%.gz}
spike_root=${spike_fasta%.fasta}
spike_root=${spike_root%.fa}
gunzip -c $spike_fasta_gz > $spike_fasta
anno_gtf=${anno_gtf_gz%.gz}
anno_root=${anno_gtf%.gtf}
gunzip -c $anno_gtf_gz > $anno_gtf

archive_file="${genome}_${anno}_$(basename ${spike_root})_rsemIndex.tgz"
if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
    archive_file="${genome}_${gender}_${anno}_$(basename ${spike_root})_rsemIndex.tgz"
fi
echo "-- Results will be: '${archive_file}'."

echo "-- Prepare reference..."
set -x
mkdir out
rsem-prepare-reference --gtf $anno_gtf ${ref_fasta},${spike_fasta} out/rsem
set +x

# Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
echo "-- Create bam header..."
set -x
refComment="@CO\tREFID:$(basename ${ref_root})"
annotationComment="@CO\tANNID:$(basename ${anno_root})"
spikeInComment="@CO\tSPIKEID:$(basename ${spike_root})"
echo -e ${refComment} > out/rsem_bamCommentLines.txt
echo -e ${annotationComment} >> out/rsem_bamCommentLines.txt
echo -e ${spikeInComment} >> out/rsem_bamCommentLines.txt
echo `cat "out/rsem_bamCommentLines.txt"`
set +x

echo "* Tar up index..."
set -x
tar -czf $archive_file out/ # * out
set +x

echo "-- The results..."
ls -l $archive_file

