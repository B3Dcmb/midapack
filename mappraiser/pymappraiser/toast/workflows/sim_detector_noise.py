import numpy as np
from astropy import units as u
import toast
import toast.ops

from sotodlib.toast import ops as so_ops
from sotodlib.toast.workflows.job import workflow_timer

import pymappraiser.toast.ops as pto


def setup_simulate_detector_noise(operators):
    """Add commandline args and operators for simulating detector noise.

    Args:
        operators (list):  The list of operators to extend.

    Returns:
        None

    """
    operators.append(pto.MySimNoise(name="my_sim_noise", enabled=False))


@workflow_timer
def simulate_detector_noise(job, otherargs, runargs, data):
    """Simulate the intrinsic detector and readout noise.

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
        job_ops.my_sim_noise.realization = otherargs.realization

    if job_ops.my_sim_noise.enabled:
        job_ops.my_sim_noise.apply(data)
