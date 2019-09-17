#!/usr/bin/env python3
"""
Script to run bam to signals (bigwigs) step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import shlex
import subprocess
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("bam_to_signals.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)

STAR_COMMAND = """STAR --runMode inputAlignmentsFromBAM \
                --inputBAMfile {input_bam} \
                --outWigType bedGraph \
                --outWigStrand {strandedness} \
                --outWigReferencesPrefix chr"""


def main(args):
    print(args)
    star_return_code = call_star(args.bamfile, args.strandedness)

    try:
        assert star_return_code == 0
    except AssertionError:
        logger.exception("Building bedGraph had a problem, most likely out of memory.")
        sys.exit(1)

    if args.strandedness == "stranded":
        call_bg_to_bw(
            "Signal.UniqueMultiple.str1.out.bg",
            args.chrom_sizes,
            args.bamroot + "_minusAll.bw",
        )
        call_bg_to_bw(
            "Signal.Unique.str1.out.bg",
            args.chrom_sizes,
            args.bamroot + "_minusUniq.bw",
        )
        call_bg_to_bw(
            "Signal.UniqueMultiple.str2.out.bg",
            args.chrom_sizes,
            args.bamroot + "_plusAll.bw",
        )
        call_bg_to_bw(
            "Signal.Unique.str2.out.bg", args.chrom_sizes, args.bamroot + "_plusUniq.bw"
        )
    else:
        call_bg_to_bw(
            "Signal.UniqueMultiple.str1.out.bg",
            args.chrom_sizes,
            args.bamroot + "_all.bw",
        )
        call_bg_to_bw(
            "Signal.Unique.str1.out.bg", args.chrom_sizes, args.bamroot + "_uniq.bw"
        )


def call_star(input_bam, strandedness):
    command = STAR_COMMAND.format(
        input_bam=input_bam, strandedness=strandedness.capitalize()
    )
    logger.info("Running STAR command %s", command)
    return_code = subprocess.call(shlex.split(command))
    return return_code


def call_bg_to_bw(input_bg, chrom_sizes, out_fn):
    # sort bedgraph
    bedgraph_cmd = "bedSort {input_bg} {output_bg}".format(
        input_bg=input_bg, output_bg=input_bg
    )
    logger.info("Sorting bedgraph: %s", bedgraph_cmd)
    subprocess.call(shlex.split(bedgraph_cmd))
    # make bigwig
    command = "bedGraphToBigWig {input_bg} {chrom_sizes} {out_fn}".format(
        input_bg=input_bg, chrom_sizes=chrom_sizes, out_fn=out_fn
    )
    logger.info("Building bigWig: %s", command)
    subprocess.call(shlex.split(command))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bamfile", type=str, help="Input bam")
    parser.add_argument(
        "--chrom_sizes", type=str, help="chromosome sizes file the input bam"
    )
    parser.add_argument("--strandedness", type=str, choices=["stranded", "unstranded"])
    parser.add_argument(
        "--bamroot",
        type=str,
        help="""
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             """,
        default="out_bam",
    )
    args = parser.parse_args()
    main(args)
