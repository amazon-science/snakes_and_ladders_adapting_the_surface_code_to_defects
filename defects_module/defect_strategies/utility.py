"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import (
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.code_library import (
    RotatedSurfaceCode,
)


def split_ancillas_by_native_type(patch: RotatedSurfaceCode) -> tuple[set[Pos], set[Pos]]:
    """Function that splits the ancillas up into X and Z types based on the stabilizers
    of the undamaged rotated surface code patch.

    Input arguments:
        `patch`: Rotated surface code patch.

    Output arguments:
        A tuple of two sets of Pos, each corresponding to all the ancilla qubits of X
        or Z type.
    """
    x_ancillas: set[Pos] = set()
    z_ancillas: set[Pos] = set()
    for ancilla in patch.ancilla_qubits:
        stabilizer_type = patch._get_stabilizer_type(ancilla.x, ancilla.y)
        if stabilizer_type == PauliT.X:
            x_ancillas.add(ancilla)
        elif stabilizer_type == PauliT.Z:
            z_ancillas.add(ancilla)
        else:
            raise ValueError(f"Unrecognized stabilizer type {stabilizer_type}.")

    return x_ancillas, z_ancillas


def stabilizers_commute(stabilizer1: Stabilizer, stabilizer2: Stabilizer) -> bool:
    """Function that checks if two stabilizers (or two gauge checks) commute.

    Input arguments:
        `stabilizer1`: First stabilizer (or gauge check). It can be a superstabilizer.
        `stabilizer2`: Second stabilizer (or gauge check). It can be a superstabilizer.

    Output arguments:
        True or False.
    """
    # They commute if of the same type
    if stabilizer1.type == stabilizer2.type:
        return True
    # otherwise check that the number of shared qubits between each stabilizer
    # and all gauges of the other is even
    commute = all(
        [
            len([qubit for qubit in set(gauge.data_qubits) & set(stabilizer2.data_qubits)]) % 2 == 0
            for gauge in SuperStabilizer.decompose(
                [
                    stabilizer1,
                ]
            )
        ]
    )
    commute &= all(
        [
            len([qubit for qubit in set(gauge.data_qubits) & set(stabilizer1.data_qubits)]) % 2 == 0
            for gauge in SuperStabilizer.decompose(
                [
                    stabilizer2,
                ]
            )
        ]
    )
    return commute


class NoStrategyFoundError(Exception):
    """Error that is raised when a data qubit or an ancilla qubit is on the boundary
    of a patch or a window. The code is then considered to be uncorrectable.
    """

    pass
