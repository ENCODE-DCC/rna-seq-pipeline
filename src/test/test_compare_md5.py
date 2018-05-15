#!/usr/bin/env python3
'''
unittests for compare_md5.py
'''

import unittest
from unittest.mock import patch
import compare_md5
from io import BytesIO


class TestFileWithMd5(unittest.TestCase):
    def setUp(self):
        self.filepath = '/this/is/a/path/to/file.txt'

    def test_get_file_with_md5(self):
        obj = compare_md5.get_file_with_md5(self.filepath)
        self.assertTrue(isinstance(obj, compare_md5.FileWithMd5))

    def test_object_gets_paths_correctly(self):
        obj = compare_md5.FileWithMd5(self.filepath)
        self.assertEqual(obj.filepath, self.filepath)
        self.assertEqual(obj.basename, 'file.txt')

    @patch('builtins.open', return_value=BytesIO(b'foo_bar'))
    def test_calculate_md5(self, mock_open):
        obj = compare_md5.FileWithMd5(self.filepath)
        self.assertEqual(obj.md5, '5c7d96a3dd7a87850a2ef34087565a6e')
