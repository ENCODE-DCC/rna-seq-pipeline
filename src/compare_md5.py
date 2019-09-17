#!/usr/bin/env python3
"""
Script for comparing md5 sums of results to a reference json.
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import hashlib
import json
import logging
import os
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("compare_md5.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


class FileWithMd5(object):
    """Class that implements a file with its md5 sum.

    Attributes:
        filepath: The path to the file.
        basename: The basename of the filepath
        __md5: string that contains the md5 sum of the file. Calculated once
        when needed the first time, and then looked up.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.basename = os.path.basename(filepath)
        self.__md5 = None

    @property
    def md5(self):
        if self.__md5 is not None:
            return self.__md5
        else:
            self.__md5 = self.calculate_md5()
            return self.__md5

    def calculate_md5(self, chunksize=4096):
        hash_md5 = hashlib.md5()
        with open(self.filepath, "rb") as f:
            # Iter is calling f.read(chunksize) until it returns
            # the sentinel b''(empty bytes)
            for chunk in iter(lambda: f.read(chunksize), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()


def get_file_with_md5(filepath):
    return FileWithMd5(filepath)


def flatten_list(input_list):
    """Flattens a nested list.

    Args:
        input_list: A (possibly) nested list.

    Returns:
        A flattened list, preserving order.
    """

    if not input_list:
        return []
    if isinstance(input_list[0], list):
        return flatten_list(input_list[0]) + flatten_list(input_list[1:])
    else:
        return input_list[:1] + flatten_list(input_list[1:])


def main(args):
    with open(args.reference_json) as f:
        reference = json.load(f)
    with open(args.metadata_json) as f:
        metadata = json.load(f)
    output = metadata["outputs"]
    # make the output structure flat
    for key in output:
        if isinstance(output[key], list):
            output[key] = flatten_list(output[key])
        else:
            output[key] = [output[key]]

    # build the list of input files to inspect
    input_files = []
    for key in args.keys_to_inspect:
        try:
            for item in output[key]:
                input_files.append(item)
        except KeyError:
            logger.exception("Key %s was not found in the output!", key)
            sys.exit(-1)

    files_to_inspect = [get_file_with_md5(file) for file in input_files]
    files_to_inspect_md5 = {file.basename: file.md5 for file in files_to_inspect}
    md5_match_by_file = dict()
    match_overall = True
    try:
        for key in files_to_inspect_md5:
            match = reference[key] == files_to_inspect_md5[key]
            md5_match_by_file[key] = match
            match_overall &= match
    except KeyError:
        logger.exception("Key %s not found", key)
        match_overall = False
        md5_match_by_file["match_overall"] = False
    else:
        md5_match_by_file["match_overall"] = match_overall
    with open(args.outfile, "w") as f:
        json.dump(md5_match_by_file, fp=f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--keys_to_inspect", nargs="+")
    parser.add_argument("--metadata_json")
    parser.add_argument("--reference_json")
    parser.add_argument("--outfile")
    args = parser.parse_args()
    main(args)
