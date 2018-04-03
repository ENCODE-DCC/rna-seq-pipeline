#!/usr/bin/env python3
"""
Unit test module for mad_qc.py
"""

import unittest
import mad_qc


class TestMadQC(unittest.TestCase):
    def setUp(self):
        self.empty_string = ''
        self.somestring = 'somestring'
        self.fn_results = 'filename.results'
        self.fn_txt = 'filename.txt'
        self.rep1_something_results = 'rep1SOMETHING.results'
        self.rep2_something_results = 'rep2SOMETHING.results'
        self.rep1_something_text = 'rep1SOMETHING.txt'
        self.rep2_something_text = 'rep2SOMETHING.txt'
        self.rep1_variable_txt = 'rep1fgfgfgfg.txt'
        self.rep1_variable_results = 'rep1zxzxzxzxzx.results'

    def test_divide_on_common_empty_empty(self):
        self.assertEqual(
            mad_qc.divide_on_common(self.empty_string, self.empty_string),
            (['', '', ''], ['', '', '']))

    def test_divide_on_common_empty_something(self):
        self.assertEqual(mad_qc.divide_on_common(self.empty_string, self.somestring), (['', '', ''], ['', '', '']) )