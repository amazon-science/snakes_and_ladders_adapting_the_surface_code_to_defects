"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

import numpy as np
import sinter

from defects_module.translate import PhysicalCircuit


def simulate_decode(circ: PhysicalCircuit, num_shots: int = 10_000) -> np.ndarray:
    """Simulates a QEC experiment, and decodes the observable outcomes.
    Returns an array of shape (num_shots, num_obs)
    """
    ### Sample measurements
    sampler = circ.stim_circuit.compile_sampler()
    msmt_data = sampler.sample(num_shots)

    ### Convert to detectors
    converter = circ.stim_circuit.compile_m2d_converter()
    dets = converter.convert(measurements=msmt_data.astype(np.bool_), separate_observables=False)

    ### Retrieve observables
    observable_data = np.array(
        [msmt_data[:, obs].sum(axis=1) % 2 for obs in circ.observables],
        dtype=np.bool_,
    ).T

    ### Build decoder
    detector_error_model = circ.stim_circuit.detector_error_model(
        allow_gauge_detectors=False,
        approximate_disjoint_errors=True,
        decompose_errors=True,
        ignore_decomposition_failures=False,
    )

    ### Decode
    corrected_observable_data = observable_data ^ sinter.predict_observables(
        dem=detector_error_model,
        dets=dets,
        decoder="pymatching",
    )

    return corrected_observable_data


def get_logical_error_rate(
    circ: PhysicalCircuit, num_shots: int = 10_000, obs_ind: int = 0
) -> float:
    """Returns the logical statistics of the specified observable index.
    Generally, we only write memory experiments with this package, and
    therefore have only a single observable.
    """
    corrected_observable_data = simulate_decode(circ, num_shots)
    return np.sum(corrected_observable_data[:, obs_ind]) / num_shots
