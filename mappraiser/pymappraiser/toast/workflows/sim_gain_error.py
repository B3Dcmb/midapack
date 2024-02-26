# Copyright (c) 2023-2023 Simons Observatory.
# Full license can be found in the top level "LICENSE" file.
"""Simulated calibration error.
"""

import numpy as np
from astropy import units as u
import toast
import toast.ops

from sotodlib.toast import ops as so_ops
from sotodlib.toast.workflows.job import workflow_timer

import pymappraiser.toast.ops as pto


def setup_simulate_calibration_error(operators):
    """Add commandline args and operators for simulating gain error.

    Args:
        operators (list):  The list of operators to extend.

    Returns:
        None

    """
    operators.append(pto.MyGainScrambler(name="my_gainscrambler", enabled=False))


@workflow_timer
def simulate_calibration_error(job, otherargs, runargs, data):
    """Simulate calibration errors.

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

    if otherargs.realization is not None:
        job_ops.my_gainscrambler.realization = otherargs.realization

    if job_ops.my_gainscrambler.enabled:
        job_ops.my_gainscrambler.apply(data)
