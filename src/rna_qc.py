#!/usr/bin/env python3
"""
Script to run additional QC step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging

import pysam
from qc_utils import QCMetric, QCMetricRecord

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("rna_qc.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


def read_dict_from_tsv(path_to_tsv):
    """Reads tab separated file with two columns into a dict.

    Args:
        path_to_tsv: filepath that contains the tsv with two columns.

    Returns:
        result: dict with key value pairs so that the first item in a row is
        the key, and the second item is  the value.

    Raises:
        AssertionError if a line in input does not have exactly two columns.
    """

    result = {}
    with open(path_to_tsv, "rt") as f:
        for line in f:
            line_split = line.split()
            try:
                assert len(line_split) == 2, "Malformed line: {}".format(line)
            except AssertionError:
                logger.exception("Malformed line %s", line)
                raise
            result.update({line_split[0]: line_split[1]})
    return result


def get_gene_type_counts(tr_to_gene_type_map, bampath):
    """Counts reads by gene type from transcriptome alignment .bam file.

    Args:
        tr_to_gene_type_map: dict that maps transcript ids to gene types
        bampath: file path to transcriptome alignment .bam

    Returns:
        reads_by_gene_type: dict with counts of reads by gene type

    Raises:
        KeyError if transcript id (read.reference_name) that is in the .bam
        is not found in the tr_to_gene_type_map. The intended use of this
        function is to use the mapping generated from the same annotation
        that was used to build the bam, so this is an indication of mismatched
        inputs.
    """
    bamfile = pysam.AlignmentFile(bampath, "rb")
    reads_by_gene_type = {key: 0 for key in set(tr_to_gene_type_map.values())}
    reads_by_gene_type.update({"transcript_id_not_found": 0})
    for read in bamfile.fetch(until_eof=True):
        if read.is_secondary or read.is_unmapped or read.is_qcfail or read.is_duplicate:
            continue
        else:
            try:
                reference_name = read.reference_name
                gene_type = tr_to_gene_type_map[reference_name]
                reads_by_gene_type[gene_type] += 1
            except KeyError:
                logger.exception(
                    "Transcript ID %s not found in mapping!", read.reference_name
                )
                reads_by_gene_type["transcript_id_not_found"] += 1
    return reads_by_gene_type


def main(args):
    qc_record = QCMetricRecord()
    logger.info(
        "Reading transcript id to gene type mapping from %s",
        args.tr_id_to_gene_type_tsv,
    )
    tr_to_gene_type_map = read_dict_from_tsv(args.tr_id_to_gene_type_tsv)
    logger.info("Calculating gene type counts for bam %s", args.input_bam)
    gene_type_counts = get_gene_type_counts(tr_to_gene_type_map, args.input_bam)
    gene_type_counts = QCMetric("gene_type_count", gene_type_counts)
    qc_record.add(gene_type_counts)
    logger.info("Writing QC output into %s", args.output_filename)
    with open(args.output_filename, "wt") as fp:
        json.dump(qc_record.to_ordered_dict(), fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_bam", type=str, help="path to transcriptome bam")
    parser.add_argument(
        "--tr_id_to_gene_type_tsv",
        type=str,
        help="path to transcript id to gene type tsv",
    )
    parser.add_argument("--output_filename", type=str)
    args = parser.parse_args()
    main(args)
