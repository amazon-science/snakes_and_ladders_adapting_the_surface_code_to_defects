"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools

import networkx

from defects_module.base import (
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.code_library import (
    RotatedSurfaceCode,
)
from defects_module.defect_strategies.boundary_deformation.utility import BoundaryDeformation
from defects_module.defect_strategies.make_windows import Window
from defects_module.defect_strategies.utility import (
    split_ancillas_by_native_type,
)

from .make_holes import Hole


def collect_data_from_holes(
    window: Window, holes: list[Hole]
) -> tuple[set[Stabilizer], list[BoundaryDeformation]]:
    """Function that collects the properties from all the holes in the patch.
    It returns the stabilizers of the patch and all the possible boundary deformations.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `holes`: A list of `Hole` objects each storing information about a specific
        hole in the patch, such as the superstabilizers and the possible boundary
        deformations that can be applied around it if it is open.

    Output arguments:
        A tuple containing the set of stabilizers in the patch and a list of
        BoundaryDeformation objects, one for each possible boundary deformation
        (due to corner placements).
    """
    # Get the stabilizers in the patch for each ancilla that is not
    # involved in superstabilizers and is not a defect
    patch = RotatedSurfaceCode(
        window.distance, window.sw_offset, window.sw_type, window.vertical_logical
    )
    # We define the set of stabilizers in the patch
    # NB: stabilizers that involve repurposed ancillas will be
    # later promoted to the set of superstabilizers
    stabilizers: set[Stabilizer] = set()
    # We define the list of possible boundary deformations
    # for each open hole
    possible_boundary_deformations: list[list[BoundaryDeformation]] = []
    # We define the sets of disabled qubits
    data_disabled: set[Pos] = set()
    ancilla_disabled: set[Pos] = set()
    # We store the ancillas by type
    x_ancillas, z_ancillas = split_ancillas_by_native_type(patch)
    # we collect all the superstabilizers and possible boundary deformations
    for hole in holes:
        # the defective and zombie qubits are defined as disabled
        data_disabled |= hole.data_qubits_defective | hole.data_qubits_zombie
        ancilla_disabled |= hole.ancilla_qubits_defective | hole.ancilla_qubits_zombie
        # we use a hack where we store the repurposed ancillas
        # as disabled ancillas because we want to promote the
        # stabilizers of the surface code with those ancillas
        # to superstabilizers (measured every other round)
        gauges = [gauge for s in hole.superstabilizers["initial"] for gauge in s.gauges]
        x_gauges = [gauge for gauge in gauges if gauge.type == PauliT.X]
        z_gauges = [gauge for gauge in gauges if gauge.type == PauliT.Z]
        for gauges, ancillas in [
            (x_gauges, z_ancillas),
            (z_gauges, x_ancillas),
        ]:
            for gauge in gauges:
                if gauge.only_ancilla in ancillas:
                    continue
                ancilla_disabled.add(gauge.only_ancilla)
        # if the hole is closed
        if hole.is_closed:
            # if the hole is closed then we add its superstabilizers
            # to the sets of X or Z superstabilizers
            for superstabilizer in hole.superstabilizers["initial"]:
                if any(superstabilizer.ancilla):
                    stabilizers.add(superstabilizer)
        # if the hole is open
        else:
            # we add the possible boundary deformations
            possible_boundary_deformations.append(hole.boundary_deformations)
    # then we add the undamaged stabilizers
    for stab in patch.stabilizers:
        # skip if the ancilla is damaged
        if stab.only_ancilla in ancilla_disabled:
            continue
        # skip if any data qubit is damaged
        if any([data_qubit in data_disabled for data_qubit in stab.data_qubits]):
            continue
        # add to the set of undamaged stabilizers
        stabilizers.add(stab)
    # We define all combinations of boundary deformations for
    # the open holes
    boundary_deformations = [
        BoundaryDeformation.merge(list(combination))
        for combination in itertools.product(*possible_boundary_deformations)
    ]
    return stabilizers, boundary_deformations


def stabilizers_connected_components(stabilizers: set[Stabilizer]) -> list[set[Pos]]:
    """Function that returns all the groups of qubits that are connected together via
    stabilizers. Each group technically forms a potential patch. We return the list
    of these groups (potential patches).

    The function works as follows:
    1) Given a set of stabilizers, we produce a graph whose vertices are the qubits.
    A gate in each stabilizer results in an edge in the graph.
    2) We returns the connected components of this graph.

    Input arguments:
        `stabilizers`: A set of all the stabilizers in the patch.

    Output arguments:
        A list of set of Pos (ancilla and data qubits) where each set corresponds to
        a group of qubits that are connected together via stabilizers.
    """
    # We create a graph with edges corresponding to the gates in the stabilizers
    stabilizer_graph = networkx.Graph()
    for stab in stabilizers:
        for ancilla in stab.ancilla:
            for data in stab.data_qubits:
                stabilizer_graph.add_edge(ancilla, data)
    # We return the list of connected components in the graph, corresponding to
    # all the groups of qubits are connected together
    return list(networkx.connected_components(stabilizer_graph))


def get_repurposed_stabilizers(
    superstabilizers: set[SuperStabilizer], undamaged_stabilizers: set[Stabilizer]
) -> set[Stabilizer]:
    """Function that returns the repurposed stabilizers to build efficient checks.
    We promote those stabilizers to superstabilizers in the algorithm since they are
    measured every other round.

    Input arguments:
        `superstabilizers`: A set of all the superstabilizers in the patch.
        `undamaged_stabilizers`: A set of all the undamaged stabilizers in the patch.

    Output arguments:
        A set of stabilizers: they are the updated stabilizers of the patch.
    """
    repurposed: dict[Pos, list[Stabilizer | None]] = dict()
    # keep track of which ancillas are already part of superstabilizers
    # but we dont't need to keep the gauges, they will not be updated
    # i.e. promoted to superstabilizers
    for gauge in SuperStabilizer.decompose([stab for stab in superstabilizers]):
        repurposed[gauge.only_ancilla] = [None]
    # now we add all the stabilizers to the dictionary
    for gauge in undamaged_stabilizers:
        if gauge.only_ancilla not in repurposed:
            repurposed[gauge.only_ancilla] = [gauge]
        else:
            repurposed[gauge.only_ancilla].append(gauge)
    # the repurposed gauges are all those where the ancilla is used twice
    return set(
        gauge for gauges in repurposed.values() if len(gauges) == 2 for gauge in gauges if gauge
    )


class NoPatchFoundError(Exception):
    """Error that is raised when a data qubit or an ancilla qubit is on the boundary
    of a patch or a window. The code is then considered to be uncorrectable.
    """

    pass


def update_patch_stabilizers(
    stabilizers: set[Stabilizer], boundary_deformation: BoundaryDeformation
) -> tuple[set[Stabilizer], set[SuperStabilizer]]:
    """Function that updates the stabilizers of the patch given a boundary
    deformation (defined as a `BoundaryDeformation` object). It returns the
    updated undamaged stabilizers and superstabilizers of the patch separately.

    The function works as follows:
    1) It finds the largest group of qubits that are connected together via
    stabilizers to form a patch. If it finds none, then it raises an error.
    2) It promotes all the stabilizers involving repurposed ancillas to
    superstabilizers since they will be measured every other round and therefore
    treated the same way as the other superstabilizers. They are however still
    technically stabilizers.

    Input arguments:
        `stabilizers`: Stabilizers of the patch obtained from the `Hole` objects.
        `boundary_deformation`: A `BoundaryDeformation` object containing all
        the information about the boundary deformation to be applied.

    Output arguments:
        A tuple containing:
        1) A set of undamaged stabilizers.
        2) A set of superstabilizers.
    """
    # First we update the set of stabilizers
    updated_stabilizers = stabilizers.copy()
    updated_stabilizers |= boundary_deformation.superstabilizers
    updated_stabilizers |= boundary_deformation.boundary_checks
    updated_superstabilizers = set(
        stab for stab in updated_stabilizers if isinstance(stab, SuperStabilizer)
    )
    updated_undamaged_stabilizers = updated_stabilizers - updated_superstabilizers
    # We remove disconnected components
    connected_components = stabilizers_connected_components(updated_stabilizers)
    if not connected_components:
        # no stabilizers are found in the patch: this is a failure
        # and therefore we skip the current strategy
        raise NoPatchFoundError(
            "Strategy dropped: Did not find connected stabilizers to form a patch."
        )
    subpatch = sorted(connected_components, key=len)[-1]
    updated_undamaged_stabilizers = set(
        stab for stab in updated_undamaged_stabilizers if stab.only_ancilla in subpatch
    )
    # We find all repurposed ancillas and update the superstabilizers
    repurposed_stabilizers = get_repurposed_stabilizers(
        updated_superstabilizers, updated_undamaged_stabilizers
    )
    updated_superstabilizers |= set(SuperStabilizer([gauge]) for gauge in repurposed_stabilizers)
    updated_undamaged_stabilizers -= repurposed_stabilizers
    return (
        updated_undamaged_stabilizers,
        updated_superstabilizers,
    )
