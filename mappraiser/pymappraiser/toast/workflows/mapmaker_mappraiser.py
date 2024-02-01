# Copyright (c) 2023-2023 Simons Observatory.
# Full license can be found in the top level "LICENSE" file.
"""Mapmaking with the MADAM destriper.
"""

import numpy as np
from astropy import units as u
import toast
import toast.ops
from toast.observation import default_values as defaults

from sotodlib.toast import ops as so_ops
from sotodlib.toast.workflows.job import workflow_timer

import pymappraiser.toast.ops as pto


def setup_mapmaker_mappraiser(parser, operators):
    """Add commandline args and operators for the MAPPRAISER mapmaker.

    Args:
        parser (ArgumentParser):  The parser to update.
        operators (list):  The list of operators to extend.

    Returns:
        None

    """
    parser.add_argument(
        "--ref",
        required=False,
        default="run0",
        help="Reference that is added to the name of the output maps.",
    )

    if pto.mappraiser.available():
        operators.append(pto.Mappraiser(name="mappraiser", enabled=True))


@workflow_timer
def mapmaker_mappraiser(job, otherargs, runargs, data):
    """Run the MAPPRAISER mapmaker.

    Args:
        job (namespace):  The configured operators and templates for this job.
        otherargs (namespace):  Other commandline arguments.
        runargs (namespace):  Job related runtime parameters.
        data (Data):  The data container.

    Returns:
        None

    """
    # Configured operators for this job
    job_ops = job.operators

    if pto.mappraiser.available() and job_ops.mappraiser.enabled:
        job_ops.mappraiser.params["path_output"] = otherargs.out_dir
        job_ops.mappraiser.params["ref"] = otherargs.ref
        job_ops.mappraiser.pixel_pointing = job.pixels_final
        job_ops.mappraiser.stokes_weights = job.weights_final
        job_ops.mappraiser.apply(data)
