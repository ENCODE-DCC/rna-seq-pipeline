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
from bisect import insort
from pysam import AlignmentFile

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
    """Container that holds the qc metric dict and the "master" key
    of the said qc metric.
    """

    def __init__(self, qc_metric_name, qc_metric_dict):
        if not isinstance(qc_metric_dict, dict):
            raise TypeError('QCMetric data must be a dict.')
        self._qc_metric_name = qc_metric_name
        self._qc_metric_dict = qc_metric_dict

    @property
    def qc_metric_dict(self):
        return self._qc_metric_dict

    @property
    def qc_metric_name(self):
        return self._qc_metric_name

    def __lt__(self, other):
        return self.qc_metric_name < other.qc_metric_name

    def __eq__(self, other):
        return self.qc_metric_name == other.qc_metric_name

    def __repr__(self):
        return 'QCMetric(%r, %r)' % (self.qc_metric_name, self.qc_metric_dict)


class QCMetricRecord(object):
    """Container that holds QCMetrics in sorted order.

    Attributes:
        metrics: list of metrics, kept sorted by the name of metrics
    """

    def __init__(self, metrics=None):
        if metrics is None:
            self._metrics = []

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
            qc_metric.qc_metric_name)
        insort(self._metrics, qc_metric)

    def __iter__(self):
        """Iterating QCMetricRecord is iterating over metrics.
        """
        return iter(self.metrics)

    def __repr__(self):
        """Like __iter__, __repr__ is delegated to metrics.
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
            assert len(line_split) == 2, 'Malformed line: {}'.format(line)
            result.update({line_split[0]: line_split[1]})
    return result
