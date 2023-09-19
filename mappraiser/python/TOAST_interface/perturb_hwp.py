import os
from time import time

import numpy as np
import traitlets
from astropy import units as u

from toast import rng
from toast.observation import default_values as defaults
from toast.timing import Timer, function_timer
from toast.traits import Int, Quantity, Unicode, Bool, trait_docs
from toast.utils import GlobalTimers, Logger, Timer, dtype_to_aligned, name_UID
from toast.ops.operator import Operator


class MyPerturbHWP(Operator):
    """Operator that adds irregularities to HWP rotation"""

    # Class traits

    API = Int(0, help="Internal interface version for this operator")

    times = Unicode(defaults.times, help="Observation shared key for timestamps")

    hwp_angle = Unicode(
        defaults.hwp_angle,
        allow_none=True,
        help="Observation shared key for HWP angle",
    )

    drift_sigma = Quantity(
        None,
        allow_none=True,
        help="1-sigma relative change in spin rate, such as 0.01 / hour",
    )

    time_sigma = Quantity(
        None,
        allow_none=True,
        help="1-sigma difference between real and nominal time stamps",
    )

    realization = Int(0, allow_none=False, help="Realization index")

    output_dir = Unicode(
        ".",
        help="Write products to this directory",
    )

    write_angle = Bool(
        False,
        help="If True, write the old and new hwp_angle buffers",
    )

    @function_timer
    def _exec(self, data, detectors=None, **kwargs):
        t0 = time()
        log = Logger.get()

        for trait in ("times", "hwp_angle"):
            if getattr(self, trait) is None:
                msg = f"You must set the '{trait}' trait before calling exec()"
                raise RuntimeError(msg)

        for iobs, obs in enumerate(data.obs):
            offset = obs.local_index_offset
            nlocal = obs.n_local_samples
            ntotal = obs.n_all_samples

            # Get an RNG seed
            key1 = self.realization * 1543343 + obs.telescope.uid
            key2 = obs.session.uid
            counter1 = 0

            # The times and hwp_angle are shared among columns of the process
            # grid.  Only the first process row needs to modify the data.
            if (
                obs.shared.comm_type(self.times) != "column"
                or obs.shared.comm_type(self.hwp_angle) != "column"
            ):
                msg = f"obs {obs.name}: expected shared fields {self.times} and "
                msg += f"{self.hwp_angle} to be on the column communicator."
                raise RuntimeError(msg)

            if obs.comm_col_rank == 0:
                times = obs.shared[self.times].data
                hwp_angle = obs.shared[self.hwp_angle].data

                # We are in the first process row.  In our RNG generation,
                # "counter2" corresponds to the sample index.  If there are
                # multiple processes in the grid row, start our RNG stream
                # at the first sample on this process.
                counter2 = obs.local_index_offset

                time_delta = times[-1] - times[0]

                # Simulate timing error (jitter)
                if self.time_sigma is None:
                    time_error = 0
                else:
                    component = 0
                    rngdata = rng.random(
                        times.size,
                        sampler="gaussian",
                        key=(key1, key2 + component),
                        counter=(counter1, counter2),
                    )
                    time_error = np.array(rngdata) * self.time_sigma.to_value(u.s)
                new_times = times + time_error
                if np.any(np.diff(new_times) <= 0):
                    raise RuntimeError("Simulated timing error causes time travel")

                # Simulate rate drift
                unwrapped = np.unwrap(hwp_angle)
                median_step = np.median(np.diff(unwrapped))
                if np.abs(median_step) < 1e-10:
                    # This was a stepped HWP, not continuously rotating
                    msg = "Don't know how to perturb a stepped HWP. "
                    msg += f"Median step size is {np.degrees(median_step)} deg"
                    raise ValueError(msg)
                nominal_rate = (unwrapped[-1] - unwrapped[0]) / time_delta
                if self.drift_sigma is None:
                    begin_rate = nominal_rate
                    accel = 0
                else:
                    # This random number is for the uniform drift across the whole
                    # observation.  All processes along the row of the grid should
                    # use the same value here.
                    counter2 = 0
                    component = 1
                    rngdata = rng.random(
                        1,
                        sampler="gaussian",
                        key=(key1, key2 + component),
                        counter=(counter1, counter2),
                    )
                    sigma = self.drift_sigma.to_value(1 / u.s) * time_delta
                    drift = rngdata[0] * sigma
                    begin_rate = nominal_rate * (1 - drift)
                    end_rate = nominal_rate * (1 + drift)
                    accel = (end_rate - begin_rate) / time_delta

                # Now calculcate the HWP angle subject to jitter and drift
                t = new_times - new_times[0]
                new_angle = 0.5 * accel * t**2 + begin_rate * t + hwp_angle[0]
                new_angle %= 2 * np.pi
                
                # Save the data
                if self.write_angle:
                    fname = os.path.join(self.output_dir, f"{iobs:04d}_{obs.uid}_hwp")
                    np.savez(
                        fname, hwp_angle=obs.shared[self.hwp_angle].data, new_angle=new_angle
                    )
            else:
                new_angle = None

            # Set the new HWP angle values
            obs.shared[self.hwp_angle].set(new_angle, offset=(0,), fromrank=0)

    def _finalize(self, data, **kwargs):
        return

    def _requires(self):
        return {
            "shared": [
                self.times,
                self.hwp_angle,
            ]
        }

    def _provides(self):
        return dict()
