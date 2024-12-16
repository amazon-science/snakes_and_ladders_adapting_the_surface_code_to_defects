"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from defects_module.base import NoiseModel, PauliT
from defects_module.code_library import RotatedSurfaceCode
from defects_module.translate import PhysicalCircuit


def standard_noise(p: float) -> NoiseModel:
    """Produces a standard noise model."""
    return {
        "R": [("DEPOLARIZE1", p)],
        "M": [("M", p)],
        "H": [("DEPOLARIZE1", p)],
        # "I": [("DEPOLARIZE1", p)],
        "CX": [("DEPOLARIZE2", p)],
    }


def memory_program(
    patch: RotatedSurfaceCode,
    noise_model: NoiseModel,
    pauli: PauliT = PauliT.X,
    rounds: int = 0,
) -> PhysicalCircuit:
    """Makes a physical circuit which prepares a logical qubit in the given
    Pauli basis, does several rounds of checks, and then measures in the
    same Pauli basis.
    """
    circ = PhysicalCircuit(noise_model)

    patch.initialize(circ, pauli)
    for _ in range(rounds):
        patch.check(circ)
    patch.measure(circ, pauli)

    return circ
