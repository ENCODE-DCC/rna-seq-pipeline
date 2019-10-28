#!/usr/bin/env python3
"""
Script to run madQC step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import shlex
import subprocess

from rna_qc import QCMetric, QCMetricRecord

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("mad_qc.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)

MADQC_CMD = "Rscript {path_to_madR} {quants_1} {quants_2}"


def remove_quantfile_extensions(quant_fn):
    first_extension_start_index = quant_fn.find(".")
    if first_extension_start_index == -1:
        return quant_fn
    else:
        return quant_fn[:first_extension_start_index]


def main(args):
    run_cmd = MADQC_CMD.format(
        path_to_madR=args.MAD_R_path, quants_1=args.quants1, quants_2=args.quants2
    )
    quant_basename1 = remove_quantfile_extensions(os.path.basename(args.quants1))
    quant_basename2 = remove_quantfile_extensions(os.path.basename(args.quants2))
    plot_output_filename = "{basename_1}-{basename_2}_mad_plot.png".format(
        basename_1=quant_basename1, basename_2=quant_basename2
    )
    # capture the output string from the run
    logger.info("Running madQC command %s", run_cmd)
    mad_output = subprocess.check_output(shlex.split(run_cmd))
    os.rename("MAplot.png", plot_output_filename)
    qc_record = QCMetricRecord()
    mad_r_metric = json.loads(mad_output.decode())
    mad_r_metric_obj = QCMetric("MAD.R", mad_r_metric)
    qc_record.add(mad_r_metric_obj)
    qc_output_fn = "{basename_1}-{basename_2}_mad_qc_metrics.json".format(
        basename_1=quant_basename1, basename_2=quant_basename2
    )
    with open(qc_output_fn, "w") as f:
        json.dump(qc_record.to_ordered_dict(), f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quants1", type=str, help="first quantification file from RSEM"
    )
    parser.add_argument(
        "--quants2", type=str, help="second quantification file from RSEM"
    )
    parser.add_argument("--MAD_R_path", type=str, help="path to MAD.R")
    args = parser.parse_args()
    main(args)
