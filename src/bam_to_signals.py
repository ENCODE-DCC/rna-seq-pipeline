#!/usr/bin/env python3
"""
Script to run bam to signals (bigwigs) step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import subprocess
import shlex

STAR_COMMAND = '''STAR --runMode inputAlignmentsFromBAM \
                --inputBAMfile {input_bam} \
                --outWigType bedGraph \
                --outWigStrand {strandedness} \
                --outFileNamePrefix ./Signal/ \
                --outWigReferencesPrefix chr'''


def main(args):
    call_star(args.bamfile, args.strandedness)
    if args.strandedness == 'stranded':
        call_bg_to_bw('Signal.UniqueMultiple.str1.out.bg', args.chrom_sizes,
                      args.bamroot + '_minusAll.bw')
        call_bg_to_bw('Signal.Unique.str1.out.bg', args.chrom_sizes,
                      args.bamroot + '_minusUniq.bw')
        call_bg_to_bw('Signal.UniqueMultiple.str2.out.bg', args.chrom_sizes,
                      args.bamroot + '_plusAll.bw')
        call_bg_to_bw('Signal.Unique.str2.out.bg', args.chrom_sizes,
                      args.bamroot + '_plusUniq.bw')
    else:
        call_bg_to_bw('Signal.UniqueMultiple.str1.out.bg', args.chrom_sizes,
                      args.bamroot + '_all.bw')
        call_bg_to_bw('Signal.Unique.str1.out.bg', args.chrom_sizes,
                      args.bamroot + '_uniq.bw')


def call_star(input_bam, strandedness):
    command = STAR_COMMAND.format(
        input_bam=input_bam, strandedness=strandedness)
    subprocess.call(shlex.split(command))


def call_bg_to_bw(input_bg, chrom_sizes, out_fn):
    command = 'bedGraphToBigWig {input_bg} {chrom_sizes}'.format(
        input_bg=input_bg, chrom_sizes=chrom_sizes)
    with open(out_fn, 'wb') as f:
        subprocess.call(shlex.split(command), stdout=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bamfile', type=str, nargs=1, help='Input bam')
    parser.add_argument(
        '--chrom_sizes',
        type=str,
        nargs=1,
        help='chromosome sizes file the input bam')
    parser.add_argument(
        '--strandedness', type=str, choices=['stranded', 'unstranded'])
    parser.add_argument(
        '--bamroot',
        type=str,
        help='''
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             ''',
        default='out_bam')
    args = parser.parse_args()
    main(args)
