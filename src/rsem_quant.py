#!/usr/bin/env python3
"""
Script to run alignment(mapping) step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

from align import make_modified_TarInfo
import argparse
import logging
import os
import re
import shlex
import subprocess
import tarfile

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler('rsem_quant.log')
filehandler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
logger.addHandler(filehandler)

RSEM_COMMAND = '''rsem-calculate-expression --bam \
--estimate-rspd \
--calc-ci \
--seed {rnd_seed} \
-p {ncpus} \
--no-bam-output \
--ci-memory {ramGB}000 \
--forward-prob {fwd_prob} \
{paired_end} \
{anno_bam} \
rsem_index/rsem \
{bam_root}_rsem'''


def strand_to_fwd_prob(strand):
    if strand == 'forward':
        return 1
    if strand == 'reverse':
        return 0
    if strand == 'unstranded':
        return 0.5
    raise ValueError('Strand must be forward, reverse or unstranded')


def format_endedness(endedness):
    if endedness == 'paired':
        return '--paired-end'
    else:
        return ''


def main(args):
    remove_bam_from_end_re = re.compile('\.bam$')
    bam_root = remove_bam_from_end_re.sub('', os.path.basename(args.anno_bam))
    with tarfile.open(args.rsem_index, 'r:gz') as archive:
        archive.extractall(
            '.', members=make_modified_TarInfo(archive, 'rsem_index'))
    rsem_call = shlex.split(
        RSEM_COMMAND.format(
            rnd_seed=args.rnd_seed,
            ncpus=args.ncpus,
            ramGB=args.ramGB,
            fwd_prob=strand_to_fwd_prob(args.read_strand),
            paired_end=format_endedness(args.endedness),
            anno_bam=args.anno_bam,
            bam_root=bam_root))
    logger.info('Running RSEM command %s', ' '.join(rsem_call))
    subprocess.call(rsem_call)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--rsem_index', type=str, help='RSEM index gzipped tar')
    parser.add_argument(
        '--anno_bam', type=str, help='STAR alignment to annotation.')
    parser.add_argument('--endedness', type=str, choices=['paired', 'single'])
    parser.add_argument(
        '--read_strand',
        type=str,
        choices=['forward', 'reverse', 'unstranded'])
    parser.add_argument(
        '--rnd_seed', type=int, help='random seed', default=12345)
    parser.add_argument('--ncpus', type=int, help='number of cpus available')
    parser.add_argument('--ramGB', type=int, help='memory available in GB')
    args = parser.parse_args()
    main(args)
