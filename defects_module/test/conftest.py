"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

import pytest

import defects_module as dm


@pytest.fixture
def memory_circuit():
    """Logical memory circuit for a distance 5 surface code."""
    patch = dm.RotatedSurfaceCode(5)
    noise_model = dm.standard_noise(0.001)
    circ = dm.memory_program(patch, noise_model, dm.PauliT.X, rounds=5)

    return circ


@pytest.fixture
def damaged_patch():
    """Distance 5 surface code with a single data defect in the center."""
    undamaged_patch = dm.RotatedSurfaceCode(5)
    q = dm.Pos(5, 5)

    x_gauges, z_gauges = set(), set()
    undamaged_stabilizers = set()
    for stab in undamaged_patch.stabilizers:
        if q in stab.data_qubits:
            if stab.type == dm.PauliT.X:
                x_gauges.add(
                    dm.Stabilizer(stab.ancilla, tuple(p for p in stab.pauli if p.qubit != q))
                )
            else:
                z_gauges.add(
                    dm.Stabilizer(stab.ancilla, tuple(p for p in stab.pauli if p.qubit != q))
                )
        else:
            undamaged_stabilizers.add(stab)

    return dm.DefectiveSurfaceCode(
        patch=undamaged_patch,
        effective_distance=None,
        undamaged_stabilizers=undamaged_stabilizers,
        x_superstabilizers=set([dm.SuperStabilizer(x_gauges)]),
        z_superstabilizers=set([dm.SuperStabilizer(z_gauges)]),
        ancilla_defects=set(),
        data_defects=set([q]),
        link_defects=set(),
    )
