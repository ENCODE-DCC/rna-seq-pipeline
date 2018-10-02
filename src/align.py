#!/usr/bin/env python3
"""
Script to run alignment(mapping) step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

from abc import ABC, abstractmethod
import argparse
import logging
import os
import shlex
import subprocess
import tarfile

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler('align.log')
filehandler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
logger.addHandler(filehandler)


def make_aligner(args):
    if args.aligner == 'star':
        if args.endedness == 'single':
            logger.info('Creating a single-ended aligner.')
            return SingleEndedStarAligner(args.fastqs, args.ncpus, args.ramGB,
                                          args.indexdir, args.bamroot)
        elif args.endedness == 'paired':
            logger.info('Creating a paired-end aligner.')
            return PairedEndStarAligner(args.fastqs, args.ncpus, args.ramGB,
                                        args.indexdir, args.bamroot)


def make_modified_TarInfo(archive, target_dir=''):
    """
    Input: archive TarFile object
    Input: target_dir string
    Returns a list of modified TarInfo objects that are files, whose
    extraction path gets modified into target_dir.
    Example: archive index.tgz contents are out (directory)
    out/file1.txt, out/file2.txt. Extracting the files
    into cwd under directory my_output (instead of ./out/)
    can be done with this:
    with tarfile.open('index.tgz', 'r:gz') as archive:
        archive.extractall('.', members=make_modified_TarInfo(archive,
        'my_output'))
    """

    members = []
    for member in archive.getmembers():
        if member.isfile():
            member.name = os.path.join(target_dir,
                                       os.path.basename(member.name))
            members.append(member)
    return members


def get_flagstats(input_path, output_path):
    command = 'samtools flagstat {infile}'.format(infile=input_path)
    logger.info('Getting flagstats for %s', input_path)
    process = subprocess.run(shlex.split(command), stdout=subprocess.PIPE)
    with open(output_path, 'w') as f:
        f.write(process.stdout.decode())


class StarAligner(ABC):
    """
    Abstract base class that gathers aspects common to both PE and SE
    Star aligning jobs.
    """

    def __init__(self, ncpus, ramGB, indexdir, bamroot):
        self.ncpus = ncpus
        self.ramGB = ramGB
        self.indexdir = indexdir
        self.bamroot = bamroot

    def run(self):
        logger.info('running STAR command %s', ' '.join(self.command))
        subprocess.call(self.command)

    @property
    @abstractmethod
    def command_string(self):
        pass

    @abstractmethod
    def format_command_string(self):
        pass

    def post_process(self):
        cwd = os.getcwd()
        os.rename(
            os.path.join(cwd, 'Aligned.sortedByCoord.out.bam'),
            os.path.join(cwd, self.bamroot + '_genome.bam'))
        os.rename(
            os.path.join(cwd, 'Log.final.out'),
            os.path.join(cwd, self.bamroot + '_Log.final.out'))
        rsem_check_cmd = 'rsem-sam-validator {bam_to_check}'.format(
            bam_to_check='Aligned.toTranscriptome.out.bam')
        rsem_output = subprocess.check_output(shlex.split(rsem_check_cmd))
        # rsem validator exits with 0 whether the check passes or not
        # for this reason we check if the output ends in 'is valid!'
        # the other possibility is 'is not valid!'
        rsem_valid = rsem_output.decode().strip().split('\n')[-1].endswith(
            'is valid!')
        if rsem_valid:
            logger.info('Transcriptome bam is already rsem-sorted.')
            os.rename(
                os.path.join(cwd, 'Aligned.toTranscriptome.out.bam'),
                os.path.join(cwd, self.bamroot + '_anno.bam'))
        else:
            logger.info('Transcriptome bam is not rsem-sorted.')
            rsem_sort_cmd = 'convert-sam-for-rsem {input} {output}'.format(
                input='Aligned.toTranscriptome.out.bam',
                output=self.bamroot + '_anno')
            logger.info('Running %s', rsem_sort_cmd)
            subprocess.call(shlex.split(rsem_sort_cmd))


class SingleEndedStarAligner(StarAligner):

    command_string = '''STAR --genomeDir {indexdir} \
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

    def __init__(self, fastqs, ncpus, ramGB, indexdir, bamroot):
        super().__init__(ncpus, ramGB, indexdir, bamroot)
        self.input_fastq = fastqs[0]
        self.command = shlex.split(
            self.format_command_string(type(self).command_string))

    def format_command_string(self, input_string):
        cmd = input_string.format(
            infastq=self.input_fastq,
            ncpus=self.ncpus,
            ramGB=self.ramGB,
            indexdir=self.indexdir)
        return cmd


class PairedEndStarAligner(StarAligner):

    command_string = '''STAR --genomeDir {indexdir} \
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

    def __init__(self, fastqs, ncpus, ramGB, indexdir, bamroot):
        super().__init__(ncpus, ramGB, indexdir, bamroot)
        self.fastq_read1 = fastqs[0]
        self.fastq_read2 = fastqs[1]
        self.command = shlex.split(
            self.format_command_string(type(self).command_string))

    def format_command_string(self, input_string):
        cmd = input_string.format(
            read1_fq_gz=self.fastq_read1,
            read2_fq_gz=self.fastq_read2,
            ncpus=self.ncpus,
            ramGB=self.ramGB,
            indexdir=self.indexdir)
        return cmd


def main(args):
    logger.info('running %s-end %s aligner', args.endedness, args.aligner)
    with tarfile.open(args.index, 'r:gz') as archive:
        archive.extractall()
    aligner = make_aligner(args)
    aligner.run()
    aligner.post_process()
    cwd = os.getcwd()
    genome_bam_path = os.path.join(cwd, args.bamroot + '_genome.bam')
    anno_bam_path = os.path.join(cwd, args.bamroot + '_anno.bam')
    genome_flagstat_path = os.path.join(cwd,
                                        args.bamroot + '_genome_flagstat.txt')
    anno_flagstat_path = os.path.join(cwd, args.bamroot + '_anno_flagstat.txt')
    get_flagstats(genome_bam_path, genome_flagstat_path)
    get_flagstats(anno_bam_path, anno_flagstat_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastqs', nargs='+', help='Input gzipped fastq(s)')
    parser.add_argument(
        '--aligner',
        type=str,
        choices=['star'],
        help='this is exists for the purpose of extension',
        required=True,
        default='star')
    parser.add_argument(
        '--index',
        type=str,
        help='Path to aligner index tar.gz archive.',
        required=True)
    parser.add_argument(
        '--indexdir',
        type=str,
        help='Directory to extract index to.',
        default='out')
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
