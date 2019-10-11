"""
unittests for rsem_quant.py
"""

import unittest
from unittest.mock import patch

import rsem_quant
from pandas import DataFrame


class TestFunctions(unittest.TestCase):
    def test_strand_to_fwd_prob_values(self):
        self.assertEqual(rsem_quant.strand_to_fwd_prob("forward"), 1)
        self.assertEqual(rsem_quant.strand_to_fwd_prob("reverse"), 0)
        self.assertEqual(rsem_quant.strand_to_fwd_prob("unstranded"), 0.5)

    def test_strand_to_fwd_prob_raises(self):
        with self.assertRaises(KeyError):
            rsem_quant.strand_to_fwd_prob("foo bar baz")

    def test_format_endedness(self):
        self.assertEqual(rsem_quant.format_endedness("paired"), "--paired-end")
        self.assertEqual(rsem_quant.format_endedness(None), "")
        self.assertEqual(rsem_quant.format_endedness("foo bar"), "")

    @patch(
        "pandas.read_csv",
        return_value=DataFrame([{"TPM": 0.1}, {"TPM": 2.0}, {"TPM": 0.8}]),
    )
    def test_calculate_number_of_genes_detected(self, mock_read):
        self.assertEqual(
            rsem_quant.calculate_number_of_genes_detected("path/to/tsv", 1), 1
        )
        self.assertEqual(
            rsem_quant.calculate_number_of_genes_detected("path/to/tsv", 0.7), 2
        )
        self.assertEqual(
            rsem_quant.calculate_number_of_genes_detected("/path/to/tsv", 3), 0
        )
