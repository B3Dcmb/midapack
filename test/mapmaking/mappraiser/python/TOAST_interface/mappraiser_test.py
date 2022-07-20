# Unit test for the toast3 / mappraiser interface
# This code is mostly taken from toast/tests/ops_madam.py

import os

import numpy as np
import numpy.testing as nt
from astropy import units as u

from . import new_mappraiser
from toast import ops as ops
from toast.noise import Noise
from toast.observation import default_values as defaults
from toast.pixels import PixelData, PixelDistribution
from toast.vis import set_matplotlib_backend
from toast.tests._helpers import (
    create_fake_sky,
    create_outdir,
    create_satellite_data,
    fake_flags,
)
from toast.tests.mpi import MPITestCase


class MappraiserTest(MPITestCase):
    def setUp(self):
        fixture_name = os.path.splitext(os.path.basename(__file__))[0]
        self.outdir = create_outdir(self.comm, fixture_name)
        np.random.seed(123456)

    def test_mappraiser_interface(self):
        if not new_mappraiser.available():
            print("libmappraiser not available, skipping tests")
            return

        # Create a fake satellite data set for testing
        data = create_satellite_data(self.comm)

        # Create some detector pointing matrices
        detpointing = ops.PointingDetectorSimple()
        pixels = ops.PixelsHealpix(
            nside=64,
            detector_pointing=detpointing,
            create_dist="pixel_dist",
        )
        pixels.apply(data)
        weights = ops.StokesWeights(
            mode="IQU",
            hwp_angle=defaults.hwp_angle,
            detector_pointing=detpointing,
        )
        weights.apply(data)

        # Create fake polarized sky pixel values locally
        create_fake_sky(data, "pixel_dist", "fake_map")

        # Scan map into timestreams
        scanner = ops.ScanMap(
            det_data=defaults.det_data,
            pixels=pixels.pixels,
            weights=weights.weights,
            map_key="fake_map",
        )
        scanner.apply(data)

        # if data.comm.world_rank == 0:
        #     set_matplotlib_backend()
        #     import matplotlib.pyplot as plt
        #
        #     ob = data.obs[0]
        #     det = ob.local_detectors[0]
        #     xdata = ob.shared["times"].data
        #     ydata = ob.detdata["signal"][det]
        #
        #     fig = plt.figure(figsize=(12, 8), dpi=72)
        #     ax = fig.add_subplot(1, 1, 1, aspect="auto")
        #     ax.plot(
        #         xdata,
        #         ydata,
        #         marker="o",
        #         c="red",
        #         label="{}, {}".format(ob.name, det),
        #     )
        #     # cur_ylim = ax.get_ylim()
        #     # ax.set_ylim([0.001 * (nse.NET(det) ** 2), 10.0 * cur_ylim[1]])
        #     ax.legend(loc=1)
        #     plt.title("Sky Signal")
        #     savefile = os.path.join(
        #         self.outdir, "signal_sky_{}_{}.pdf".format(ob.name, det)
        #     )
        #     plt.savefig(savefile)
        #     plt.close()

        # Create an uncorrelated noise model from focalplane detector properties
        default_model = ops.DefaultNoiseModel(noise_model="noise_model")
        default_model.apply(data)

        # Simulate noise from this model
        sim_noise = ops.SimNoise(noise_model="noise_model", out="noise")
        sim_noise.apply(data)

        # Make fake flags
        # fake_flags(data)

        # if data.comm.world_rank == 0:
        #     set_matplotlib_backend()
        #     import matplotlib.pyplot as plt
        #
        #     ob = data.obs[0]
        #     det = ob.local_detectors[0]
        #     xdata = ob.shared["times"].data
        #     ydata = ob.detdata["signal"][det]
        #
        #     fig = plt.figure(figsize=(12, 8), dpi=72)
        #     ax = fig.add_subplot(1, 1, 1, aspect="auto")
        #     ax.plot(
        #         xdata,
        #         ydata,
        #         marker="o",
        #         c="red",
        #         label="{}, {}".format(ob.name, det),
        #     )
        #     # cur_ylim = ax.get_ylim()
        #     # ax.set_ylim([0.001 * (nse.NET(det) ** 2), 10.0 * cur_ylim[1]])
        #     ax.legend(loc=1)
        #     plt.title("Sky + Noise Signal")
        #     savefile = os.path.join(
        #         self.outdir, "signal_sky-noise_{}_{}.pdf".format(ob.name, det)
        #     )
        #     plt.savefig(savefile)
        #     plt.close()

        # Compute timestream rms

        # rms = dict()
        # for ob in data.obs:
        #     rms[ob.name] = dict()
        #     for det in ob.local_detectors:
        #         flags = np.array(ob.shared[defaults.shared_flags])
        #         flags |= ob.detdata[defaults.det_flags][det]
        #         good = flags == 0
        #         # Add an offset to the data
        #         ob.detdata[defaults.det_data][det] += 500.0
        #         rms[ob.name][det] = np.std(ob.detdata[defaults.det_data][det][good])

        # if data.comm.world_rank == 0:
        #     set_matplotlib_backend()
        #     import matplotlib.pyplot as plt
        #
        #     ob = data.obs[0]
        #     det = ob.local_detectors[0]
        #     xdata = ob.shared["times"].data
        #     ydata = ob.detdata["signal"][det]
        #
        #     fig = plt.figure(figsize=(12, 8), dpi=72)
        #     ax = fig.add_subplot(1, 1, 1, aspect="auto")
        #     ax.plot(
        #         xdata,
        #         ydata,
        #         marker="o",
        #         c="red",
        #         label="{}, {}".format(ob.name, det),
        #     )
        #     # cur_ylim = ax.get_ylim()
        #     # ax.set_ylim([0.001 * (nse.NET(det) ** 2), 10.0 * cur_ylim[1]])
        #     ax.legend(loc=1)
        #     plt.title("Sky + Noise + Offset Signal")
        #     savefile = os.path.join(
        #         self.outdir, "signal_sky-noise-offset_{}_{}.pdf".format(ob.name, det)
        #     )
        #     plt.savefig(savefile)
        #     plt.close()

        # Run mappraiser on this

        # Mappraiser assumes constant sample rate- just get it from the noise model for
        # the first detector.
        sample_rate = data.obs[0]["noise_model"].rate(data.obs[0].local_detectors[0])

        pars = {}
        pars["fsample"] = sample_rate
        pars["nside"] = pixels.nside
        pars["path_output"] = self.outdir
        pars["solver"] = 0
        pars["precond"] = 0
        pars["Z_2lvl"] = 0
        pars["pointing_commflag"] = 6
        pars["tol"] = 1e-6
        pars["maxiter"] = 500
        pars["enlFac"] = 1
        pars["ortho_alg"] = 1
        pars["bs_red"] = 0
        pars["Lambda"] = 4096
        pars["ref"] = "run0"

        # FIXME: add a view here once our test data includes it

        mappraiser = new_mappraiser.Mappraiser(
            params=pars,
            det_data=defaults.det_data,
            pixel_pointing=pixels,
            stokes_weights=weights,
            noise_model="noise_model",
            noise_name="noise",
            copy_groups=2,
            purge_det_data=False,
            restore_det_data=False,
        )
        mappraiser.apply(data)

        # if data.comm.world_rank == 0:
        #     set_matplotlib_backend()
        #     import matplotlib.pyplot as plt
        #
        #     ob = data.obs[0]
        #     det = ob.local_detectors[0]
        #     xdata = ob.shared["times"].data
        #     ydata = ob.detdata["destriped"][det]
        #
        #     fig = plt.figure(figsize=(12, 8), dpi=72)
        #     ax = fig.add_subplot(1, 1, 1, aspect="auto")
        #     ax.plot(
        #         xdata,
        #         ydata,
        #         marker="o",
        #         c="red",
        #         label="{}, {}".format(ob.name, det),
        #     )
        #     # cur_ylim = ax.get_ylim()
        #     # ax.set_ylim([0.001 * (nse.NET(det) ** 2), 10.0 * cur_ylim[1]])
        #     ax.legend(loc=1)
        #     plt.title("Destriped Signal")
        #     savefile = os.path.join(
        #         self.outdir, "signal_destriped_{}_{}.pdf".format(ob.name, det)
        #     )
        #     plt.savefig(savefile)
        #     plt.close()

        # for ob in data.obs:
        #     for det in ob.local_detectors:
        #         flags = np.array(ob.shared[defaults.shared_flags])
        #         flags |= ob.detdata[defaults.det_flags][det]
        #         good = flags == 0
        #         check_rms = np.std(ob.detdata["destriped"][det][good])
        #         # print(f"check_rms = {check_rms}, det rms = {rms[ob.name][det]}")
        #         self.assertTrue(0.9 * check_rms < rms[ob.name][det])

        del data
        return
