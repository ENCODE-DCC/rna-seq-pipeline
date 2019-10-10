#!/usr/bin/env python3
"""
Script to run kallisto quantification step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import os
import shlex
import subprocess
import sys
from abc import ABC, abstractmethod

from align import concatenate_files

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("kallisto_quant.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)


class KallistoQuant(ABC):
    """Abstract base class for kallisto quant objects.

    Attributes:
        command_template: formatting template for the run command as a class
            attribute. Override in subclasses
        output_dir: Is the directory where output gets written relative to the
            current working directory.
        number_of_threads: Int that determines number of threads used.
        strandedness: String: forward, reverse or unstranded
    """

    def __init__(
        self, path_to_index, output_dir, number_of_threads, strandedness, out_prefix
    ):
        self.path_to_index = path_to_index
        self.output_dir = output_dir
        self.number_of_threads = number_of_threads
        self.strandedness_direction = self.parse_strandedness(strandedness)
        self.out_prefix = out_prefix

    @property
    @abstractmethod
    def command_template(self):
        """class constant: command template for kallisto quant
        """
        pass

    def format_command_template(self, **kwargs):
        """Method that adds concrete values to the command template
        """
        command = self.command_template.format(**kwargs)
        return command

    @property
    def command(self):
        return self._command

    def run(self):
        logger.info("Running kallisto command: %s", " ".join(self.command))
        subprocess.call(self.command)

    @staticmethod
    def parse_strandedness(strandedness):
        """Transform the strandedness from input to kallisto format.

        Args:
            strandedness: string that is either forward, reverse or unstranded

        Returns:
            string
            forward -> --fr-stranded
            reverse -> --rf-stranded
            unstranded -> ""

        Raises:
            KeyError if input is not forward, reverse or unstranded
        """
        strandedness_dict = {
            "forward": "--fr-stranded",
            "reverse": "--rf-stranded",
            "unstranded": "",
        }

        return strandedness_dict[strandedness]

    def rename_output(self):
        kallisto_out_absolute = os.path.join(os.getcwd(), self.output_dir)
        os.rename(
            os.path.join(kallisto_out_absolute, "abundance.tsv"),
            os.path.join(kallisto_out_absolute, self.out_prefix + "_abundance.tsv"),
        )


class KallistoQuantSingleEnd(KallistoQuant):
    """Class that sets up kallisto quant run in single-ended mode:

    Attributes:
        fragment_length: Int that determines the fragment length of the library
        sd_of_fragment_length: Double that determines the standard deviation of
             the fragment_length.
        fastq: list with the path to input fastq.
    """

    command_template = """kallisto quant -i {path_to_index} \
    -o {output_dir} \
    -t {number_of_threads} \
    --plaintext \
    {strandedness_direction} \
    --single \
    -l {fragment_length} \
    -s {sd_of_fragment_length} \
    {fastq}"""

    def __init__(
        self,
        path_to_index,
        output_dir,
        number_of_threads,
        strandedness,
        fragment_length,
        sd_of_fragment_length,
        fastqs,
        out_prefix,
    ):
        super().__init__(
            path_to_index, output_dir, number_of_threads, strandedness, out_prefix
        )
        self.fragment_length = fragment_length
        self.sd_of_fragment_length = sd_of_fragment_length
        # we should get only one input fastq
        try:
            assert len(fastqs) == 1
            self.fastq = fastqs[0]
        except AssertionError:
            logger.exception("More than one input fastqs in single-ended mode.")
            sys.exit(1)
        self._command = shlex.split(
            self.format_command_template(
                path_to_index=self.path_to_index,
                output_dir=self.output_dir,
                number_of_threads=self.number_of_threads,
                strandedness_direction=self.strandedness_direction,
                fragment_length=self.fragment_length,
                sd_of_fragment_length=self.sd_of_fragment_length,
                fastq=self.fastq,
            )
        )


class KallistoQuantPairedEnd(KallistoQuant):
    """Class that sets up running kallisto quant in paired-end mode

    Attributes:
        see superclass
    """

    command_template = """ kallisto quant -i {path_to_index} \
    -o {output_dir} \
    -t {number_of_threads} \
    --plaintext \
    {strandedness_direction} \
    {fastq1} {fastq2}"""

    def __init__(
        self,
        path_to_index,
        output_dir,
        number_of_threads,
        strandedness,
        fastqs,
        out_prefix,
    ):
        super().__init__(
            path_to_index, output_dir, number_of_threads, strandedness, out_prefix
        )
        self.fastq1 = fastqs[0]
        self.fastq2 = fastqs[1]
        self._command = shlex.split(
            self.format_command_template(
                path_to_index=self.path_to_index,
                output_dir=self.output_dir,
                number_of_threads=self.number_of_threads,
                strandedness_direction=self.strandedness_direction,
                fastq1=self.fastq1,
                fastq2=self.fastq2,
            )
        )


def main(args):
    merged_R2 = None
    if len(args.fastqs_R1) > 1:
        merged_R1 = concatenate_files(args.fastqs_R1)
    else:
        merged_R1 = args.fastqs_R1[0]
    if args.endedness == "paired" and len(args.fastqs_R2) > 1:
        merged_R2 = concatenate_files(args.fastqs_R2)
    elif args.endedness == "paired" and len(args.fastqs_R2) == 1:
        merged_R2 = args.fastqs_R2[0]
    fastqs = [merged_R1]

    if merged_R2 and args.endedness == "paired":
        fastqs.append(merged_R2)
    if args.endedness == "paired":
        kallisto_quantifier = KallistoQuantPairedEnd(
            args.path_to_index,
            args.output_dir,
            args.number_of_threads,
            args.strandedness,
            fastqs,
            args.out_prefix,
        )
    if args.endedness == "single":
        kallisto_quantifier = KallistoQuantSingleEnd(
            args.path_to_index,
            args.output_dir,
            args.number_of_threads,
            args.strandedness,
            args.fragment_length,
            args.sd_of_fragment_length,
            fastqs,
            args.out_prefix,
        )
    kallisto_quantifier.run()
    kallisto_quantifier.rename_output()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fastqs_R1", nargs="+", help="Input gzipped fastq(s) belonging to read1"
    )
    parser.add_argument(
        "--fastqs_R2", nargs="*", help="Input gzipped fastq(s) belonging to read2"
    )
    parser.add_argument(
        "--number_of_threads", type=int, help="Number of parallel threads to use."
    )
    parser.add_argument(
        "--strandedness", type=str, choices=["forward", "reverse", "unstranded"]
    )
    parser.add_argument("--path_to_index", type=str, help="Path to kallisto index.")
    parser.add_argument(
        "--output_dir", type=str, help="Output directory path", default="kallisto_out"
    )
    parser.add_argument(
        "--endedness",
        type=str,
        choices=["single", "paired"],
        help="Endedness of the experiment. Choices are single or paired",
    )
    parser.add_argument(
        "--fragment_length",
        type=int,
        help="""In single-ended mode fragment length of the library. \
                Not needed in paired end mode where it can be inferred \
                from the data.""",
    )
    parser.add_argument(
        "--sd_of_fragment_length",
        type=float,
        help="In single-ended mode, the standard deviation of fragment length",
    )
    parser.add_argument(
        "--out_prefix",
        type=str,
        help="Prefix for output. Prefix foo outputs to foo_abundances.tsv",
        default="",
    )
    args = parser.parse_args()
    main(args)
