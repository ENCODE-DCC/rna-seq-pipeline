import unittest

import merge_annotation


class TestMergeAnnotation(unittest.TestCase):
    """Test formatting annotation, tRNA and spikein files
       for the purpose of merging them in correct format"""

    def setUp(self):
        self.good_tRNA_line = 'chr1\tENSEMBL\ttRNAscan\t7930279\t7930348\t.\t-\t.\tgene_id "33323"; transcript_id "33323"; gene_type "tRNAscan"; gene_status "NULL"; gene_name "33323"; transcript_type "tRNAscan"; transcript_status "NULL"; transcript_name "33323"; level 3;\n'
        self.empty_string = ""
        self.bad_tRNA_line = 'chr1ENSEMBLtRNAscan79302797930348.-.gene_id "33323"; transcript_id "33323"; gene_type "tRNAscan";'
        self.unchanged_string = "this line will stay the same"
        self.fasta_string_with_stuff_in_left = (
            " 1\tgg\n >NAME\nATTC\n>NAME2\nCACACA\n\n\n"
        )
        self.fasta_string = ">NAME\nATTC\n>NAME2\nCACACA\n\n\n"

    def test_replace_on_good_line(self):
        self.assertEqual(
            merge_annotation.replace_nth_position_with(
                self.good_tRNA_line, 2, "exon", "\t"
            ),
            'chr1\tENSEMBL\texon\t7930279\t7930348\t.\t-\t.\tgene_id "33323"; transcript_id "33323"; gene_type "tRNAscan"; gene_status "NULL"; gene_name "33323"; transcript_type "tRNAscan"; transcript_status "NULL"; transcript_name "33323"; level 3;\n',
        )

    def test_replace_on_empty_line(self):
        self.assertRaises(
            IndexError,
            merge_annotation.replace_nth_position_with,
            self.empty_string,
            2,
            "replacement",
            " ",
        )

    def test_replace_on_bad_line(self):
        self.assertRaises(
            IndexError,
            merge_annotation.replace_nth_position_with,
            self.bad_tRNA_line,
            2,
            "replacement",
            "\t",
        )

    def test_replace_nochange(self):
        self.assertEqual(
            self.unchanged_string,
            merge_annotation.replace_nth_position_with(
                self.unchanged_string, 0, "this", " "
            ),
        )

    def test_strip_left_until_and_including_stripping(self):
        self.assertEqual(
            merge_annotation.strip_left_until_and_including(
                self.fasta_string_with_stuff_in_left, ">"
            ),
            "NAME\nATTC\n>NAME2\nCACACA\n\n\n",
        )

    def test_strip_left_until_and_including_sentinel_not_found(self):
        self.assertRaises(
            ValueError, merge_annotation.strip_left_until_and_including, "aaa", "b"
        )

    def test_get_fasta_tokens_empty(self):
        self.assertRaises(
            ValueError, merge_annotation.get_fasta_tokens, self.empty_string
        )

    def test_get_fasta_tokens_notfasta(self):
        self.assertRaises(
            ValueError, merge_annotation.get_fasta_tokens, self.unchanged_string
        )

    def test_get_fasta_tokens_good_input(self):
        self.assertEqual(
            merge_annotation.get_fasta_tokens(self.fasta_string),
            ["NAME\nATTC\n", "NAME2\nCACACA\n\n\n"],
        )

    def test_remove_whitespace_empty(self):
        self.assertEqual(
            merge_annotation.remove_whitespace(self.empty_string), self.empty_string
        )

    def test_remove_whitespace_no_whitespace(self):
        self.assertEqual(merge_annotation.remove_whitespace("aaa"), "aaa")

    def test_remove_whitespace_with_whitespace(self):
        self.assertEqual(merge_annotation.remove_whitespace("\taaa aaa\n "), "aaaaaa")


if __name__ == "__main__":
    unittest.main()
