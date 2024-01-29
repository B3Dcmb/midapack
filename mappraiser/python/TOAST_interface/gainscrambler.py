# Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import re

import traitlets

from toast import rng
from toast.observation import default_values as defaults
from toast.timing import function_timer
from toast.traits import Bool, Float, Int, Unicode, List, trait_docs
from toast.utils import Logger
from toast.ops.operator import Operator


@trait_docs
class MyGainScrambler(Operator):
    """Apply random gain errors to detector data.

    This operator draws random gain errors from a given distribution and
    applies them to the specified detectors.
    """

    # Class traits

    API = Int(0, help="Internal interface version for this operator")

    det_data_names = List(
        trait=Unicode,
        default_value=[defaults.det_data, "noise"],
        help="Observation detdata key to apply the gain error to",
    )

    pattern = Unicode(
        f".*",
        allow_none=True,
        help="Regex pattern to match against detector names. Only detectors that "
        "match the pattern are scrambled.",
    )
    center = Float(1, allow_none=False, help="Gain distribution center")

    sigma = Float(1e-3, allow_none=False, help="Gain distribution width")

    realization = Int(0, allow_none=False, help="Realization index")

    component = Int(0, allow_none=False, help="Component index for this simulation")

    coherent_gain_error = Bool(
        False,
        allow_none=False,
        help="Apply same gain error (of `sigma` percent) to all detectors. Use `pattern` "
        "to select e.g. the B detector of the detector pairs.",
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @function_timer
    def _exec(self, data, detectors=None, **kwargs):
        log = Logger.get()

        if self.pattern is None:
            pat = None
        else:
            pat = re.compile(self.pattern)

        for obs in data.obs:
            # Get the detectors we are using for this observation
            dets = obs.select_local_detectors(detectors)
            if len(dets) == 0:
                # Nothing to do for this observation
                continue

            comm = obs.comm.comm_group
            rank = obs.comm.group_rank

            sindx = obs.session.uid
            telescope = obs.telescope.uid

            focalplane = obs.telescope.focalplane

            # key1 = realization * 2^32 + telescope * 2^16 + component
            key1 = self.realization * 4294967296 + telescope * 65536 + self.component
            key2 = sindx
            counter1 = 0
            counter2 = 0

            dets_present_list = [
                set(obs.detdata[name].detectors) for name in self.det_data_names
            ]

            for det in dets:
                # Test the detector pattern
                if pat is not None and pat.match(det) is None:
                    if rank == 0:
                        msg = f"Skipping detector '{det}'"
                        log.debug(msg)
                    continue

                detindx = focalplane[det]["uid"]
                counter1 = detindx

                if self.coherent_gain_error:
                    gain = 1 + self.sigma
                else:
                    rngdata = rng.random(
                        1,
                        sampler="gaussian",
                        key=(key1, key2),
                        counter=(counter1, counter2),
                    )

                    gain = self.center + rngdata[0] * self.sigma

                for name, dets_present in zip(self.det_data_names, dets_present_list):
                    if det in dets_present:
                        obs.detdata[name][det] *= gain

        return

    def _finalize(self, data, **kwargs):
        return

    def _requires(self):
        req = {
            "meta": list(),
            "shared": list(),
            "detdata": [self.det_data],
            "intervals": list(),
        }
        return req

    def _provides(self):
        prov = {
            "meta": list(),
            "shared": list(),
            "detdata": list(),
        }
        return prov
