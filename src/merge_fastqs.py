"""
Script for merging fastqs in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__license__ = "MIT"

import argparse
import logging
import pathlib
import shutil


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("merge_fastqs.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s | %(levelname)s | %(name)s: %(message)s')
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


def main(args):
    # the output merged fastq will be simply named as the first of the input fastqs prefixed with "merged_"
    first_fastq = pathlib.Path(args.fastqs[0])
    # write into the current working directory
    merged_path = pathlib.Path("merged_" + first_fastq.name)
    logger.info("starting to write into %s" % merged_path)
    with open(merged_path, "wb") as out:
        for file in args.fastqs:
            with open(file, "rb") as add_on:
                logger.info("adding %s" % file)
                shutil.copyfileobj(add_on, out)
                logger.info("file %s successfully added" % file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastqs", nargs="+", help="fastqs to merge")
    args = parser.parse_args()
    main(args)
