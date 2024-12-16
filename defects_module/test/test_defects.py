"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

import defects_module as dm


def get_error_model(patch: dm.RotatedSurfaceCode, error_rate: float, pauli: dm.PauliT, rounds: int):
    """Function to build the QEC program for a memory experiment.
    Returns the Stim detector error model.
    """
    circ = dm.memory_program(patch, dm.standard_noise(error_rate), pauli, rounds)
    return circ.stim_circuit.detector_error_model()


def check_patch(
    d: int,
    data_defects: set[dm.Pos],
    ancilla_defects: set[dm.Pos],
    link_defects: set[tuple[dm.Pos, dm.Pos]],
    window_lims: tuple[int, int, int, int],
    effective_distance: tuple[int, int],
    superstabilizers: set[dm.Stabilizer],
):
    """Builds code patch, searches windows, adjusts stabilizers/superstabilizers.
    Then, checks against expected.
    """
    undamaged_patch = dm.RotatedSurfaceCode(d)
    _, damaged_patches = dm.DefectiveSurfaceCode.make_defective_patches(
        undamaged_patch, ancilla_defects, data_defects, link_defects, first_super_type=dm.PauliT.X
    )
    damaged_patch = damaged_patches[0]
    ancillas_in_undamaged_stabilizers = set(
        [stabilizer.only_ancilla for stabilizer in damaged_patch.undamaged_stabilizers]
    )
    ancillas_in_superstabilizers = set(
        [
            ancilla
            for superstabilizer in damaged_patch.first_super | damaged_patch.second_super
            for ancilla in superstabilizer.ancilla
        ]
    )
    assert ancillas_in_superstabilizers.isdisjoint(ancillas_in_undamaged_stabilizers)
    assert damaged_patch.effective_distance == effective_distance
    assert damaged_patch.first_super | damaged_patch.second_super == superstabilizers
    assert damaged_patch.window_lims == tuple(window_lims)
    dz, dx = undamaged_patch.distance
    get_error_model(undamaged_patch, 1e-3, dm.PauliT.X, rounds=2 * dx)
    get_error_model(undamaged_patch, 1e-3, dm.PauliT.Z, rounds=2 * dz)


def test_single_data_defect():
    d = 3
    ancilla_defects = set([])
    data_defects = set([dm.Pos(3, 3)])
    link_defects = set([])
    window_lims = (0, 6, 0, 6)
    effective_distance = (2, 2)
    superstabilizers = set(
        [
            dm.SuperStabilizer(
                [
                    dm.Stabilizer(
                        dm.Pos(2, 4),
                        [
                            dm.Pauli(dm.Pos(1, 3), dm.PauliT.X),
                            dm.Pauli(dm.Pos(1, 5), dm.PauliT.X),
                            dm.Pauli(dm.Pos(3, 5), dm.PauliT.X),
                        ],
                    ),
                    dm.Stabilizer(
                        dm.Pos(4, 2),
                        [
                            dm.Pauli(dm.Pos(3, 1), dm.PauliT.X),
                            dm.Pauli(dm.Pos(5, 1), dm.PauliT.X),
                            dm.Pauli(dm.Pos(5, 3), dm.PauliT.X),
                        ],
                    ),
                ]
            ),
            dm.SuperStabilizer(
                [
                    dm.Stabilizer(
                        dm.Pos(2, 2),
                        [
                            dm.Pauli(dm.Pos(1, 1), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(1, 3), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(3, 1), dm.PauliT.Z),
                        ],
                    ),
                    dm.Stabilizer(
                        dm.Pos(4, 4),
                        [
                            dm.Pauli(dm.Pos(3, 5), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(5, 3), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(5, 5), dm.PauliT.Z),
                        ],
                    ),
                ]
            ),
        ]
    )
    check_patch(
        d,
        data_defects,
        ancilla_defects,
        link_defects,
        window_lims,
        effective_distance,
        superstabilizers,
    )


def test_single_ancilla_defect():
    d = 3
    ancilla_defects = set([dm.Pos(2, 4)])
    data_defects = set([])
    link_defects = set([])
    window_lims = (0, 6, 0, 6)
    effective_distance = (3, 3)
    superstabilizers = set(
        [
            dm.base.SuperStabilizer(
                [
                    dm.Stabilizer(
                        dm.Pos(2, 2),
                        [dm.Pauli(dm.Pos(1, 3), dm.PauliT.X), dm.Pauli(dm.Pos(3, 3), dm.PauliT.X)],
                    ),
                    dm.Stabilizer(
                        dm.Pos(2, 6),
                        [dm.Pauli(dm.Pos(1, 5), dm.PauliT.X), dm.Pauli(dm.Pos(3, 5), dm.PauliT.X)],
                    ),
                ]
            ),
            dm.base.SuperStabilizer(
                [
                    dm.Stabilizer(
                        dm.Pos(0, 4),
                        [dm.Pauli(dm.Pos(1, 3), dm.PauliT.Z), dm.Pauli(dm.Pos(1, 5), dm.PauliT.Z)],
                    ),
                    dm.Stabilizer(
                        dm.Pos(4, 4),
                        [
                            dm.Pauli(dm.Pos(3, 3), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(3, 5), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(5, 3), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(5, 5), dm.PauliT.Z),
                        ],
                    ),
                ]
            ),
            dm.base.SuperStabilizer(
                [
                    dm.Stabilizer(
                        dm.Pos(2, 2),
                        [
                            dm.Pauli(dm.Pos(1, 1), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(3, 1), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(1, 3), dm.PauliT.Z),
                            dm.Pauli(dm.Pos(3, 3), dm.PauliT.Z),
                        ],
                    )
                ]
            ),
        ]
    )
    check_patch(
        d,
        data_defects,
        ancilla_defects,
        link_defects,
        window_lims,
        effective_distance,
        superstabilizers,
    )


def assert_patch_is_valid(patch: dm.DefectiveSurfaceCode):
    """Checks that a defective surface code patch is valid."""
    assert patch.stabilizers_commute()
    assert patch.num_encoded_qubits() == 1
    assert not patch.has_frozen_qubits()
    assert patch.is_connected()


def check_all_defective_patches(
    undamaged_patch: dm.RotatedSurfaceCode,
    ancilla_defects: set[dm.Pos],
    data_defects: set[dm.Pos],
    link_defects: set[tuple[dm.Pos, dm.Pos]],
) -> list[dm.DefectiveSurfaceCode]:
    """Runs the defect handling strategy and checks that all output patches are
    valid. Returns the output patches for further case-specific checks.
    """
    filter_distance = dm.FilterDistance(
        filter_func=dm.DistanceFilterFunctions.return_all_strategies
    )
    _, damaged_patches = dm.DefectiveSurfaceCode.make_defective_patches(
        undamaged_patch,
        ancilla_defects,
        data_defects,
        link_defects,
        repurpose_ancillas=True,
        add_padding=True,
        filter_distance=filter_distance,
        first_super_type=dm.PauliT.X,
    )
    for patch in damaged_patches:
        assert_patch_is_valid(patch)
    return damaged_patches


def check_best_defective_patches(
    undamaged_patch: dm.RotatedSurfaceCode,
    ancilla_defects: set[dm.Pos],
    data_defects: set[dm.Pos],
    link_defects: set[tuple[dm.Pos, dm.Pos]],
) -> list[dm.DefectiveSurfaceCode]:
    """Runs the defect handling strategy and checks that all output patches are
    valid. Returns the output patches for further case-specific checks.
    """
    filter_distance = dm.FilterDistance(
        filter_func=dm.DistanceFilterFunctions.return_best_strategies
    )
    _, damaged_patches = dm.DefectiveSurfaceCode.make_defective_patches(
        undamaged_patch,
        ancilla_defects,
        data_defects,
        link_defects,
        repurpose_ancillas=True,
        add_padding=True,
        filter_distance=filter_distance,
        first_super_type=dm.PauliT.X,
    )
    for patch in damaged_patches:
        assert_patch_is_valid(patch)
    return damaged_patches


def run_memory_experiment(
    patch: dm.DefectiveSurfaceCode, error_rate: float = 5e-4, shots: int = 100_000
):
    """Runs X- and Z-memory experiments with the defective patch and a low
    physical error rate. Checks that the logical error rate is better than
    the physical error rate.
    """
    for pauli in (dm.PauliT.X, dm.PauliT.Z):
        circ = dm.memory_program(patch, dm.standard_noise(error_rate), pauli, rounds=2)
        p_L = dm.get_logical_error_rate(circ, shots)
        if p_L > error_rate:
            return False
    return True


def defect_properties(damaged_patch: dm.DefectiveSurfaceCode):
    """Returns the list of properties we want to test in our snapshots."""
    return (damaged_patch._stabilizers, damaged_patch.effective_distance)


def test_edge_defect(snapshot):
    """Defective data qubit along the edge of the patch."""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(5),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 5)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert run_memory_experiment(damaged_patch)


def test_corner_hole_1(snapshot):
    """Defective data qubit near the corner of the patch."""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(5),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 1)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 2
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot
    assert run_memory_experiment(damaged_patches[0])


def test_corner_hole_2(snapshot):
    """Defective data qubit in the corner of the patch.
    Here a data defect at (1, 3) creates a corner hole
    (it makes the qubit (1, 1) frozen, and forces its removal).
    """
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(5),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 3)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 3
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot
    assert run_memory_experiment(damaged_patches[0])


def test_closed_holes_touch(snapshot):
    """Interaction of two closed holes in the bulk."""
    damaged_patches = check_best_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(4, 8)]),
        data_defects=set([dm.Pos(5, 5)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert run_memory_experiment(damaged_patch)


def test_edge_hole_touches_closed_hole_1(snapshot):
    """Edge hole interactions with an Auger hole. Case 1: they merge"""
    damaged_patches = check_best_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 9), dm.Pos(3, 7)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert (
        len([stab for stab in damaged_patch.stabilizers if isinstance(stab, dm.SuperStabilizer)])
        == 0
    )


def test_edge_hole_touches_closed_hole_2(snapshot):
    """Edge hole interactions with an Auger hole. Case 2: they remain separate"""
    damaged_patches = check_best_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 7), dm.Pos(3, 5)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert (
        len([stab for stab in damaged_patch.stabilizers if isinstance(stab, dm.SuperStabilizer)])
        == 2
    )


def test_edge_hole_touches_closed_hole_3(snapshot):
    """Edge hole interactions with an Auger hole. Case 3: they remain separate."""
    damaged_patches = check_best_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 9), dm.Pos(3, 5)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert (
        len([stab for stab in damaged_patch.stabilizers if isinstance(stab, dm.SuperStabilizer)])
        == 2
    )


def test_edge_hole_touches_closed_hole_4(snapshot):
    """Edge hole interactions with two Auger holes. One merges, one remains separate"""
    damaged_patches = check_best_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([]),
        data_defects=set([dm.Pos(1, 7), dm.Pos(3, 5), dm.Pos(3, 9)]),
        link_defects=set([]),
    )
    assert len(damaged_patches) == 1
    damaged_patch = damaged_patches[0]
    assert defect_properties(damaged_patch) == snapshot
    assert (
        len([stab for stab in damaged_patch.stabilizers if isinstance(stab, dm.SuperStabilizer)])
        == 2
    )


def test_edge_hole_touches_closed_hole_5(snapshot):
    """Edge hole interactions with an Auger hole. Case 5: Edge hole has an ancilla defect"""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(2, 6)]),
        data_defects=set([dm.Pos(1, 5), dm.Pos(5, 7)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_edge_hole_touches_closed_hole_6(snapshot):
    """Edge hole interactions with an Auger hole. Case 6: Closed hole has an ancilla defect"""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(2, 4)]),
        data_defects=set([dm.Pos(5, 7)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_edge_hole_touches_closed_hole_7(snapshot):
    """Edge hole interactions with an Auger hole. Case 7: Closed hole has an ancilla defect"""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(4, 8)]),
        data_defects=set([dm.Pos(5, 7)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_edge_hole_touches_closed_hole_8(snapshot):
    """Edge hole interactions with an Auger hole. Case 8: Edge hole has an ancilla defect"""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(2, 8)]),
        data_defects=set([dm.Pos(1, 9), dm.Pos(5, 7)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_edge_hole_touches_closed_hole_9(snapshot):
    """Edge hole interactions with an Auger hole. Case 9"""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(2, 6)]),
        data_defects=set([dm.Pos(1, 5)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_corner_hole_touches_closed_hole(snapshot):
    """Cases involving a corner open hole."""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(5),
        ancilla_defects=set([dm.Pos(4, 4)]),
        data_defects=set([dm.Pos(1, 1), dm.Pos(1, 3)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_edge_hole_becomes_corner_hole(snapshot):
    """An edge hole becomes a corner hole after cleaning gauges."""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(9),
        ancilla_defects=set([dm.Pos(2, 8)]),
        data_defects=set([dm.Pos(5, 11), dm.Pos(7, 13), dm.Pos(9, 15)]),
        link_defects=set([]),
    )
    assert [defect_properties(damaged_patch) for damaged_patch in damaged_patches] == snapshot


def test_deep_corner_hole(snapshot):
    """Corner hole cuts deeper into the patch."""
    damaged_patches = check_all_defective_patches(
        undamaged_patch=dm.RotatedSurfaceCode(7),
        ancilla_defects=set([dm.Pos(8, 8)]),
        data_defects=set([dm.Pos(5, 9), dm.Pos(9, 11), dm.Pos(11, 9), dm.Pos(13, 9)]),
        link_defects=set([]),
    )
    # Make snapshot into the one with maximal distance.
    # This is because there are too many strategies with
    # the same effective distance and even if we sort them
    # we are not guaranteed of getting the patches in
    # exactly the same order every single time.
    best_patch = max(
        damaged_patches,
        key=lambda damaged_patch: (
            max(damaged_patch.effective_distance),
            sum(damaged_patch.effective_distance),
        ),
    )
    assert defect_properties(best_patch) == snapshot
