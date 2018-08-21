#!/usr/bin/env python3
"""
Script to run kallisto quantification step in ENCODE rna-seq-pipeline
"""

__author__ = 'Otto Jolanki'
__version__ = '0.1.0'
__license__ = 'MIT'

import argparse
import subprocess
import os
import shlex
from abc import ABC, abstractmethod


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

    def __init__(self, path_to_index, output_dir, number_of_threads,
                 strandedness):
        self.path_to_index = path_to_index
        self.output_dir = output_dir
        self.number_of_threads = number_of_threads
        self.strandedness_direction = self.parse_strandedness(strandedness)

    @property
    @abstractmethod
    def command_template(self):
        """class constant: command template for kallisto quant
        """
        pass

    @abstractmethod
    def format_command_template(self):
        """Method that adds concrete values to the command template
        """
        pass

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
        strandedness_dict = {'forward': '--fr-stranded',
                             'reverse': '--rf-stranded',
                             'unstranded': ''}

        return strandedness_dict[strandedness]


class KallistoQuantSingleEnd(KallistoQuant):
    """Class that runs kallisto quant in single-ended mode:

    Attributes:
        fragment_length: Int that determines the fragment length of the library
        sd_of_fragment_length: Double that determines the standard deviation of
             the fragment_length.
        fastq: list with the path to input fastq. List because parsing in wdl
            is a pain..
    """

    command_template = '''kallisto quant -i {path_to_index} \
    -o {output_dir} \
    -t {number_of_threads} \
    --plaintext \
    {strandedness_direction} \
    --single \
    -l {fragment_length} \
    -s {sd_of_fragment_length} \
    {fastq}'''

    def __init__(self, path_to_index, output_dir, number_of_threads,
                 strandedness, fragment_length, sd_of_fragment_length,
                 fastq):
        super().__init__(path_to_index, output_dir, number_of_threads,
                         strandedness)
        self.fragment_length = fragment_length
        self.sd_of_fragment_length = sd_of_fragment_length
        # we actually get only one input
        self._fastq = fastq[0]
        self.command = shlex.split(
            self.format_command_template(type(self).command_template))

    def format_command_template(self, input_template):
        command = input_template.format(path_to_index=self.path_to_index,
                                        output_dir=self.output_dir,
                                        number_of_threads=self.number_of_threads,
                                        strandedness_direction=self.strandedness_direction,
                                        fragment_length=self.fragment_length,
                                        sd_of_fragment_length=self.sd_of_fragment_length,
                                        fastq=self._fastq)
        return command
