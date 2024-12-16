"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

import defects_module as dm


def test_memory_circuit():
    """Executes a (decoded) memory experiment with non-defective patches of
    several possible distances, and checks that physical error rate 1e-3 lies
    below threshold.
    """
    distances = (
        3,
        5,
        7,
    )
    error_rate = 5e-4
    noise_model = dm.standard_noise(error_rate)
    num_shots = 100_000

    noise_model = dm.standard_noise(error_rate)
    logical_error_rates = {}
    for d in distances:
        patch = dm.RotatedSurfaceCode(d)
        circ = dm.memory_program(patch, noise_model, dm.PauliT.X, rounds=d)
        logical_error_rates[d] = dm.get_logical_error_rate(circ, num_shots)

    assert logical_error_rates[7] <= logical_error_rates[5]
    assert logical_error_rates[5] <= logical_error_rates[3]
    assert logical_error_rates[3] <= error_rate


def test_memory_circuit_defective(damaged_patch):
    """Executes a (decoded) memory experiment with a defective patch and checks
    that physical error rate 1e-3 lies below pseudo-threshold.
    """
    error_rate = 5e-4
    noise_model = dm.standard_noise(error_rate)
    num_shots = 100_000

    noise_model = dm.standard_noise(error_rate)
    circ = dm.memory_program(damaged_patch, noise_model, dm.PauliT.X, rounds=5)
    logical_error_rate = dm.get_logical_error_rate(circ, num_shots)
    assert logical_error_rate < error_rate
