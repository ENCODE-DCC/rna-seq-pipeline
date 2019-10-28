"""
unittests for rna_qc.py
"""

import unittest
from collections import OrderedDict
from io import StringIO
from unittest.mock import patch

import rna_qc


class TestQCMetric(unittest.TestCase):
    def test_type_check(self):
        with self.assertRaises(TypeError):
            rna_qc.QCMetric("name", 1)

    def test_get_name(self):
        qc_obj = rna_qc.QCMetric("a", {})
        self.assertEqual(qc_obj.name, "a")

    def test_get_content(self):
        qc_obj = rna_qc.QCMetric("_", {2: "a", 1: "b"})
        self.assertEqual(qc_obj.content, OrderedDict([(1, "b"), (2, "a")]))

    def test_less_than(self):
        smaller_obj = rna_qc.QCMetric(1, {})
        bigger_obj = rna_qc.QCMetric(2, {})
        self.assertTrue(smaller_obj < bigger_obj)

    def test_equals(self):
        first_obj = rna_qc.QCMetric("a", {})
        second_obj = rna_qc.QCMetric("a", {"x": "y"})
        self.assertTrue(first_obj == second_obj)

    def test_repr(self):
        obj = rna_qc.QCMetric("a", {1: "x"})
        self.assertEqual(obj.__repr__(), "QCMetric('a', OrderedDict([(1, 'x')]))")


class TestQCMetricRecord(unittest.TestCase):
    def setUp(self):
        self.obj_a1 = rna_qc.QCMetric("a", {1: 2})
        self.obj_a2 = rna_qc.QCMetric("a", {2: 3})
        self.obj_b = rna_qc.QCMetric("b", {3: 4})
        self.qc_record = rna_qc.QCMetricRecord()

    def test_init_from_list_not_unique(self):
        metrics = [self.obj_a1, self.obj_a2]
        with self.assertRaises(AssertionError):
            rna_qc.QCMetricRecord(metrics)

    def test_init_from_list_success(self):
        metrics = [self.obj_a1, self.obj_b]
        record = rna_qc.QCMetricRecord(metrics)
        self.assertEqual(record.metrics[0], self.obj_a1)
        self.assertEqual(record.metrics[1], self.obj_b)

    def test_add(self):
        self.assertEqual(len(self.qc_record), 0)
        self.qc_record.add(self.obj_a1)
        self.assertEqual(len(self.qc_record), 1)

    def test_add_raises_error_when_add_same_twice(self):
        self.qc_record.add(self.obj_a1)
        with self.assertRaises(AssertionError):
            self.qc_record.add(self.obj_a1)

    def test_add_raises_error_when_add_with_same_name(self):
        self.qc_record.add(self.obj_a1)
        with self.assertRaises(AssertionError):
            self.qc_record.add(self.obj_a2)

    def test_to_ordered_dict(self):
        self.qc_record.add(self.obj_a1)
        self.qc_record.add(self.obj_b)
        qc_dict = self.qc_record.to_ordered_dict()
        self.assertEqual(qc_dict, OrderedDict([("a", {1: 2}), ("b", {3: 4})]))


class TestRegularFunctions(unittest.TestCase):
    @patch("builtins.open", return_value=StringIO("file\tcontains\tbad\ttsv\n"))
    def test_read_dict_from_tsv_malformed(self, mock_open):
        with self.assertRaises(AssertionError):
            rna_qc.read_dict_from_tsv("bad.tsv")

    @patch("builtins.open", return_value=StringIO("file\tcontains\ngood\ttsv\n"))
    def test_read_dict_from_tsv_good_input(self, mock_open):
        result_dict = rna_qc.read_dict_from_tsv("good.tsv")
        self.assertEqual(result_dict["file"], "contains")
        self.assertEqual(result_dict["good"], "tsv")
