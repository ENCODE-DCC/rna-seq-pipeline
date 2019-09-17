#!/usr/bin/env python3
"""
Create GTF file from Gencode annotations (.gtf file),
Gencode tRNAs (.gtf file) and spikeins (.fasta file).
Inputs are taken in gzipped format.
Output is given in gzipped format.
"""

__author__ = "Otto Jolanki"
__version__ = "0.1"
__license__ = "MIT"

import argparse
import gzip


def replace_nth_position_with(input_string, position, replacement, separator):
    result = input_string.split(separator)
    result[position] = replacement
    return separator.join(result)


def strip_left_until_and_including(input_string, sentinel=">"):
    """
    input: string
    output: input string - everything on the left, and
            including the sentinel character
    If the sentinel character is not found in the string
    a ValueError is raised.
    """
    return input_string[input_string.index(sentinel) + 1 :]


def get_fasta_tokens(fastastring):
    """
    input: .fasta file as string
    output: list of strings, each string representing a token in .fasta
            with everything before, and including first occurence of '>'
            removed from each.
    """
    lstripped_fastastring = strip_left_until_and_including(fastastring)
    return lstripped_fastastring.split(">")


def remove_whitespace(input_string):
    return "".join(input_string.split())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # formatting template for converting the spikeins into correct .gtf format
    SPIKEIN_FORMAT = '{spikein_name}\tspikein\texon\t1\t{sequence_length}\t.\t+\t.\tgene_id "gSpikein_{spikein_name}"; transcript_id "tSpikein_{spikein_name}";\n'
    # Annotation filename
    parser.add_argument(
        "--annotation",
        action="store",
        required=True,
        help="Gencode annotation .gtf file in gzipped format.",
    )

    # tRNA filename
    parser.add_argument(
        "--tRNA",
        action="store",
        required=True,
        help="Gencode tRNA .gtf file in gzipped format.",
    )

    # spikeins filename
    parser.add_argument(
        "--spikeins",
        action="store",
        required=True,
        help="Spikein .fasta file in gzipped format.",
    )

    parser.add_argument(
        "--output_filename",
        action="store",
        required=True,
        help="Output filename, including .gz ending.",
    )

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()

    # Process the files starts here
    annotation_file = gzip.open(args.annotation, "rt")
    tRNA_file = gzip.open(args.tRNA, "rt")
    spikein_file = gzip.open(args.spikeins, "rt")

    # Only thing that needs to be done to the annotation_file is
    # removing commentlines that start with #
    # Use generators for tRNAs and annotation to save memory

    annotation_linegenerator = (
        line for line in annotation_file if not line.startswith("#")
    )
    # remove possible comment lines from the beginning of the tRNA file
    # replace tRNAscan column with exon in the result
    tRNA_nocomments = (line for line in tRNA_file if not line.startswith("#"))
    tRNA_linegenerator = (
        replace_nth_position_with(line, 2, "exon", "\t") for line in tRNA_nocomments
    )

    with spikein_file as f:
        spikein_string = f.read()

    spikein_tokens = get_fasta_tokens(spikein_string)
    # split the tokens into sequence of [name, sequence] -pairs
    spikein_pairs = (token.split("\n", 1) for token in spikein_tokens)
    # get rid of whitespace in nucleotide code, calculate length
    spikein_name_length_pairs_no_whitespace = (
        [token[0], len("".join(token[1].split()))] for token in spikein_pairs
    )
    # build a sequence of formatted spikein lines
    spikein_lines = (
        SPIKEIN_FORMAT.format(spikein_name=token[0], sequence_length=token[1])
        for token in spikein_name_length_pairs_no_whitespace
    )

    with gzip.open(args.output_filename, "wt") as f_out:
        for line in annotation_linegenerator:
            f_out.write(line)
        for line in tRNA_linegenerator:
            f_out.write(line)
        for line in spikein_lines:
            f_out.write(line)
    annotation_file.close()
    tRNA_file.close()
