#!/usr/bin/env python3
"""
Script to run alignment(mapping) step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import tarfile


def main(args):
    ''' if (args.endedness == 'paired' and len(args.fastqs) != 2):
        raise AssertionError('paired end must have 2 fastqs')
    if (args.endedness == 'single' and len(args.fastqs) != 1):
        raise AssertionError('single end must have 1 fastq')
    '''
    print('printing args:')
    print(args)


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