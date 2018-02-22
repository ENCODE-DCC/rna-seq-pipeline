#!/bin/bash -e

if [ $# -lt 6 ] || [ $# -gt 7 ]; then
    echo "usage v1: lrna_index_tophat.sh <ref_fasta_gz> <spike_in_fasta_gz> <annotation_gtf_gz> <tiny_fq_gz> <annotation_version> <genome> [<gender>]"
    echo "Creates TopHat index of reference genome, annotation and spike-ins for long-RNA-seq. Is independent of DX and encodeD."
    exit -1; 
fi
ref_fasta_gz=$1      # Reference genome assembly in gzipped fasta format.
spike_fasta_gz=$2 # All spike-ins in single gzipped fasta format.
anno_gtf_gz=$3       # Gene annotation in gzipped gtf format
tiny_fq_gz=$4
anno=$5              # Annotation (e.g. 'v24')
genome=$6            # Genome (e.g. 'GRCh38')
if [ $# -eq 7 ]; then
    gender=$7        # Gender. Values: 'female', 'male', 'XX', 'XY' will be included in names.  Otherwise, gender neutral.  
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
tiny_fq=${tiny_fq_gz%.gz}
gunzip $tiny_fq_gz

archive_file="${genome}_${anno}_${spike_root}_tophatIndex.tgz"
if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
    archive_file="${genome}_${gender}_${anno}_${spike_root}_tophatIndex.tgz"
fi
echo "-- Results will be: '${archive_file}'."

echo "-- bowtie build..."
set -x
mkdir out
#bowtie2-build --offrate 3 -f ${ref} out/$geno_prefix  ### Definitely has an effect on results
bowtie2-build -f ${ref_fasta},$spike_fasta out/$genome
# make sure the combined fa file is preserved in the archive, so that it isn't rebuilt each time
bowtie2-inspect out/$genome > out/$genome.fa
set +x

# Attempt to make bamCommentLines.txt. NOTE tabs handled by assignment.
echo "-- Create bam header..."
set -x
refComment="@CO\tREFID:${ref_root}"
annotationComment="@CO\tANNID:${anno_root}"
echo -e ${refComment} > out/${genome}_bamCommentLines.txt
echo -e ${annotationComment} >> out/${genome}_bamCommentLines.txt
if [ -n "$spike_in" ]; then
    spikeInComment="@CO\tSPIKEID:${spike_fn}"
    echo -e ${spikeInComment} >> out/${genome}_bamCommentLines.txt
fi

echo `cat "out/${genome}_bamCommentLines.txt"`
set +x

echo "-- Run a 'quicky' tophat to generate index..."
set -x
tophat --no-discordant --no-mixed -p 8 -z0 --min-intron-length 20 --max-intron-length 1000000 \
       --read-mismatches 4 --read-edit-dist 4 --max-multihits 20 --library-type fr-firststrand \
       --GTF "$anno_root".gtf --no-coverage-search \
       --transcriptome-index=out/${anno} out/${genome} tiny.fq
set +x

echo "-- Tar up results..."
set -x
tar -czvf ${archive_file} out/${genome}* out/${anno}*
set +x

echo "-- The results..."
ls -l $archive_file

