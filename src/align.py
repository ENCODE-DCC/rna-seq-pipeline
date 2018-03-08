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


class Aligner(ABC):
    def __init__(self):
        self.command = self.format_command_string(type(self).command_string)

    def run(self):
        print('running command:')
        print(self.command)
        # subprocess.Popen(self.command)

    def set_ram(self):
        pass

    @property
    def command_string(self):
        raise NotImplementedError

    @abstractmethod
    def format_command_string(self):
        pass


class SingleEndedStarAligner(Aligner):

    command_string = '''STAR --genomeDir out --readFilesIn $reads_fq_gz \
    --readFilesCommand zcat --runThreadN $ncpus --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile COfile.txt \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate  \
    --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000'''

    def __init__(self):
        super().__init__()

    def format_command_string(self, input_string):
        return input_string


def make_aligner(aligner_program, endedness):
    pass


def make_post_processor(aligner_program, endedness):
    pass


def main(args):
    ''' if (args.endedness == 'paired' and len(args.fastqs) != 2):
        raise AssertionError('paired end must have 2 fastqs')
    if (args.endedness == 'single' and len(args.fastqs) != 1):
        raise AssertionError('single end must have 1 fastq')
    '''
    '''
    aligner = make_aligner(
        aligner_program=args.aligner, endedness=args.endedness)
    aligner.set_inputs(args.fastqs)
    aligner.set_index_path(args.index)
    aligner.set_library_id(args.libraryid)
    aligner.set_bam_root(args.bamroot)
    aligner.set_ncpus(args.ncpus)
    aligner.set_ram(args.ramGB)
    aligner.run()
    post_processor = make_post_processor(
        aligner_program=args.aligner, endedness=args.endedness)
    post_processor.run()
    '''

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
        help='Library identifier which will be added to bam header.')
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