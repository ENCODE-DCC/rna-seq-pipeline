"""
unittests for compare_md5.py
"""

import unittest
from io import BytesIO
from unittest.mock import patch

import compare_md5


class TestFileWithMd5(unittest.TestCase):
    def setUp(self):
        self.filepath = "/this/is/a/path/to/file.txt"

    def test_get_file_with_md5(self):
        obj = compare_md5.get_file_with_md5(self.filepath)
        self.assertTrue(isinstance(obj, compare_md5.FileWithMd5))

    def test_object_gets_paths_correctly(self):
        obj = compare_md5.FileWithMd5(self.filepath)
        self.assertEqual(obj.filepath, self.filepath)
        self.assertEqual(obj.basename, "file.txt")

    @patch("builtins.open", return_value=BytesIO(b"foo_bar"))
    def test_calculate_md5(self, mock_open):
        obj = compare_md5.FileWithMd5(self.filepath)
        self.assertEqual(obj.md5, "5c7d96a3dd7a87850a2ef34087565a6e")


class TestRegularFunctions(unittest.TestCase):
    def test_flatten_on_empty(self):
        self.assertEqual(compare_md5.flatten_list([]), [])

    def test_flatten_on_simple_list_of_len_1(self):
        self.assertEqual(compare_md5.flatten_list([1]), [1])

    def test_flatten_on_simple_list_of_len_2(self):
        self.assertEqual(compare_md5.flatten_list([1, "a"]), [1, "a"])

    def test_flatten_on_level_one_nested_nonempty(self):
        self.assertEqual(
            compare_md5.flatten_list([[1, "a"], [2, "b"]]), [1, "a", 2, "b"]
        )

    def test_flatten_on_level_one_nested_empty(self):
        self.assertEqual(compare_md5.flatten_list([[1, "a"], []]), [1, "a"])

    def test_flatten_deeper(self):
        self.assertEqual(
            compare_md5.flatten_list([1, [2, [3, 4, [5, 6]]], 7, [8, 9]]),
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
        )

    def test_raises_TypeError_when_called_with_int(self):
        with self.assertRaises(TypeError):
            compare_md5.flatten_list(1)

    def test_raises_TypeError_when_called_with_str(self):
        with self.assertRaises(TypeError):
            compare_md5.flatten_list("foobar")
