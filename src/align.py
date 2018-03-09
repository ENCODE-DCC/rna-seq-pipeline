#!/usr/bin/env python3
"""
Script to run alignment(mapping) step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import tarfile
from abc import ABC, abstractmethod
import subprocess
import shlex
import os


def make_aligner(args):
    if args.aligner == 'star':
        if args.endedness == 'single':
            return SingleEndedStarAligner(args.fastqs, args.ncpus, args.ramGB)
        elif args.endedness == 'paired':
            return PairedEndStarAligner(args.fastqs, args.ncpus, args.ramGB)


def make_modified_TarInfo(archive, target_dir=''):
    '''
    Input: archive Tarfile in read mode
    Input: target_dir path
    Returns a list of modified TarInfo objects that are files, whose
    extraction path gets modified into target_dir.
    '''
    members = []
    for member in archive.getmembers():
        if member.isfile():
            member.name = os.path.join(target_dir,
                                       os.path.basename(member.name))
            members.append(member)
    return members


class StarAligner(ABC):
    def __init__(self):
        pass

    def run(self):
        print('running command:')
        print(self.command)
        subprocess.call(self.command)

    @property
    @abstractmethod
    def command_string(self):
        pass

    @abstractmethod
    def format_command_string(self):
        pass

    @abstractmethod
    def post_process(self):
        pass


class SingleEndedStarAligner(StarAligner):

    command_string = '''STAR --genomeDir out \
    --readFilesIn {infastq} \
    --readFilesCommand zcat \
    --runThreadN {ncpus} \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile COfile.txt \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate  \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM {ramGB}000000000'''

    def __init__(self, fastqs, ncpus, ramGB):
        self.input_fastq = fastqs[0]
        self.ncpus = ncpus
        self.ramGB = ramGB
        self.command = shlex.split(
            self.format_command_string(type(self).command_string))

    def format_command_string(self, input_string):
        cmd = input_string.format(
            infastq=self.input_fastq, ncpus=self.ncpus, ramGB=self.ramGB)
        return cmd

    def post_process(self):
        pass


class PairedEndStarAligner(StarAligner):

    command_string = '''STAR --genomeDir out \
    --readFilesIn {read1_fq_gz} {read2_fq_gz} \
    --readFilesCommand zcat \
    --runThreadN {ncpus} \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile COfile.txt \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM {ramGB}000000000'''

    def __init__(self, fastqs, ncpus, ramGB):
        self.fastq_read1 = fastqs[0]
        self.fastq_read2 = fastqs[1]
        self.ncpus = ncpus
        self.ramGB = ramGB
        self.command = shlex.split(
            self.format_command_string(type(self).command_string))

    def format_command_string(self, input_string):
        cmd = input_string.format(
            read1_fq_gz=self.fastq_read1,
            read2_fq_gz=self.fastq_read2,
            ncpus=self.ncpus,
            ramGB=self.ramGB)
        return cmd

    def post_process(self):
        pass


def main(args):
    print('running {}-end {} aligner'.format(args.endedness, args.aligner))
    with tarfile.open(args.index, 'r:gz') as archive:
        archive.extractall()
    aligner = make_aligner(args)
    aligner.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter, add_help=True)
    parser.add_argument('--fastqs', nargs='+', help='Input gzipped fastq(s)')
    parser.add_argument(
        '--aligner',
        type=str,
        choices=['star', 'tophat'],
        help='star or tophat',
        required=True)
    parser.add_argument(
        '--index',
        type=str,
        help='Path to aligner index archive.',
        required=True)
    parser.add_argument(
        '--endedness',
        type=str,
        choices=['paired', 'single'],
        help='paired or single',
        required=True)
    parser.add_argument(
        '--libraryid',
        type=str,
        help='Library identifier which will be added to bam header.',
        default='libraryID')
    parser.add_argument(
        '--bamroot',
        type=str,
        help='''
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             ''',
        default='out_bam')
    parser.add_argument(
        '--ncpus', type=int, help='Number of cpus available.', default=4)
    parser.add_argument(
        '--ramGB', type=int, help='Amount of RAM available in GB.', default=8)

    args = parser.parse_args()
    main(args)
