"""Mapmaking with the MAPPRAISER framework."""

from pymappraiser.toast import mappraiser
from sotodlib.toast.workflows.job import workflow_timer


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

    if mappraiser.available():
        operators.append(mappraiser.Mappraiser(name="mappraiser", enabled=False))


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

    if mappraiser.available() and job_ops.mappraiser.enabled:
        job_ops.mappraiser.params["path_output"] = otherargs.out_dir
        job_ops.mappraiser.params["ref"] = otherargs.ref
        job_ops.mappraiser.pixel_pointing = job.pixels_final
        job_ops.mappraiser.stokes_weights = job.weights_final
        job_ops.mappraiser.az_name = job_ops.sim_ground.azimuth
        job_ops.mappraiser.hwpangle_name = job_ops.sim_ground.hwp_angle
        job_ops.mappraiser.apply(data)
