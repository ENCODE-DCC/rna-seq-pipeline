#!/usr/bin/env python3
"""
Script to run additional QC step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import json
import logging
import pysam
from bisect import insort
from collections import OrderedDict

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler('rna_qc.log')
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


class QCMetric(object):
    """Container that holds the qc metric as OrderedDict (sorted by keys of
    the input dict) and the "master" key (name) of the said qc metric. Can be
    instantiated from a regular dict.
    """

    def __init__(self, qc_metric_name, qc_metric_dict):
        if not isinstance(qc_metric_dict, dict):
            raise TypeError('QCMetric data must be a dict.')
        self._name = qc_metric_name
        self._content = OrderedDict(
            sorted(qc_metric_dict.items(), key=lambda x: x[0]))

    @property
    def content(self):
        return self._content

    @property
    def name(self):
        return self._name

    def __lt__(self, other):
        return self.name < other.name

    def __eq__(self, other):
        return self.name == other.name

    def __repr__(self):
        return 'QCMetric(%r, %r)' % (self.name, self.content)


class QCMetricRecord(object):
    """Container that holds QCMetrics in sorted order.

    Attributes:
        metrics: list of metrics, kept sorted by the name of metrics
    """

    def __init__(self, metrics=None):
        if metrics is None:
            self._metrics = []
        else:
            # names must be unique
            names = [metric.name for metric in metrics]
            assert len(names) == len(
                set(names)), 'Names of metrics have to be unique'
            metrics.sort()
            self._metrics = metrics

    @property
    def metrics(self):
        return self._metrics

    def add(self, qc_metric):
        """Adds qc metric to the metrics, keeping it sorted by name.

        Args:
            qc_metric: QCMetric

        Returns: None

        Raises: AssertionError if a metric with same name is already in record
        """

        assert qc_metric not in self._metrics, 'Metric with name {} already in record'.format(
            qc_metric.name)
        insort(self._metrics, qc_metric)

    def to_ordered_dict(self):
        """Returns an OrderedDict with the contents.

        Returns: Ordered dict with structure as follows:
            - Ordered as the metrics is
            - Contents, assuming metrics = [qc1, qc2, qc3]:
            {
                qc1.name : qc1.content,
                qc2.name : qc2.content,
                qc3.name : qc3.content
            }
        """
        result = OrderedDict()
        for metric in self.metrics:
            result.update({metric.name: metric.content})
        return result

    def __len__(self):
        """
        Delegated to metrics.
        """
        return len(self.metrics)

    def __iter__(self):
        """
        Iterating QCMetricRecord is iterating over metrics.
        """
        return iter(self.metrics)

    def __repr__(self):
        """
        Like __iter__, __repr__ is delegated to metrics.
        """
        return self.metrics.__repr__()


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
    with open(path_to_tsv, 'rt') as f:
        for line in f:
            line_split = line.split()
            try:
                assert len(line_split) == 2, 'Malformed line: {}'.format(line)
            except AssertionError:
                logger.exception('Malformed line %s', line)
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
    bamfile = pysam.AlignmentFile(bampath, 'rb')
    reads_by_gene_type = {key: 0 for key in set(tr_to_gene_type_map.values())}
    for read in bamfile.fetch(until_eof=True):
        if read.is_secondary or read.is_unmapped or read.is_qcfail or read.is_duplicate:
            continue
        else:
            try:
                reference_name = read.reference_name
                gene_type = tr_to_gene_type_map[reference_name]
                reads_by_gene_type[gene_type] += 1
            except KeyError:
                logger.exception('Transcript ID %s not found in mapping!',
                                 read.reference_name)
    return reads_by_gene_type


def main(args):
    qc_record = QCMetricRecord()
    logger.info('Reading transcript id to gene type mapping from %s',
                args.tr_id_to_gene_type_tsv)
    tr_to_gene_type_map = read_dict_from_tsv(args.tr_id_to_gene_type_tsv)
    logger.info('Calculating gene type counts for bam %s', args.input_bam)
    gene_type_counts = get_gene_type_counts(tr_to_gene_type_map,
                                            args.input_bam)
    gene_type_counts = QCMetric('gene_type_count', gene_type_counts)
    qc_record.add(gene_type_counts)
    logger.info('Writing QC output into %s', args.output_filename)
    with open(args.output_filename, 'wt') as fp:
        json.dump(qc_record.to_ordered_dict(), fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_bam', type=str, help='path to transcriptome bam')
    parser.add_argument(
        '--tr_id_to_gene_type_tsv',
        type=str,
        help='path to transcript id to gene type tsv')
    parser.add_argument('--output_filename', type=str)
    args = parser.parse_args()
    main(args)
