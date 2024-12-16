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
from defects_module.defect_strategies.boundary_deformation.utility import (
    BoundaryDeformation,
    BoundaryDeformationString,
    BoundaryDeformationStrings,
)
from defects_module.defect_strategies.make_windows import Window

from .make_cycles import Cycle


def get_boundary_position(cycles: list[Cycle], window: Window, pauli_to_clean: PauliT) -> int:
    """Helper function that returns the position of the open hole in the patch:
        0: Left boundary
        1: Bottom boundary
        2: Right boundary
        3: Top boundary

    Input arguments:
        `cycles`: A list of `Cycle` objects associated with all the cycles
        (or contours) within the open hole to clean. (There can be more than one.)
        `window`: A `Window` object containing information about the initial window.
        `pauli_to_clean`: Pauli type of the gauge checks to clean along the boundary.

    Output arguments:
        An integer corresponding to the boundary on which the open hole is.
    """
    # define pos = 0 (left boundary), 1 (bottom boundary),
    # 2 (right boundary) and 3 (top boundary)
    qubits_on_contours = {q for cycle in cycles for q in cycle.qubits}
    if window.vertical_logical == pauli_to_clean:
        return 2 * any([node.x == window.xlims[1] - 1 for node in qubits_on_contours])
    return 1 + 2 * any([node.y == window.ylims[1] - 1 for node in qubits_on_contours])


def find_open_cycles(cycles: list[Cycle], pauli_to_clean: PauliT) -> set[int]:
    """Helper function that finds which cycles are to be considered open within
    the open hole. The rest of the cycles would be closed.

    The function works as follows:
    1) We create a graph whose vertices are the indices of the cycles within
    the open hole. The edges are determined by the pauli type of the gauge
    checks two cycles share together. If the pauli type is equal to the pauli
    type to be cleaned `pauli_to_clean` then they are considered connected.
    This is because when we clean one of them, the other will automatically
    become open as well and will have to be cleaned as well.
    2) For each pair of cycles, we check if they share gauge checks of the
    pauli type to clean. If so, we add an edge between them in the graph.
    3) For each connected component of the graph, if all the cycles in the
    component are closed, then we move to the next. If any of the cycles is
    open, then all of the cycles in the component are considered open.
    4) We return all the open cycles within the hole.

    Input arguments:
        `cycles`: A list of `Cycle` objects associated with all the cycles
        (or contours) within the open hole to clean. (There can be more than one.)
        `pauli_to_clean`: Pauli type of the gauge checks to clean along the boundary.

    Output arguments:
        A set of integers, where each integer correspong to the index of a
        cycle in the open hole that is considered open.
    """
    # create a graph of connected cycles
    graph = networkx.Graph()
    # number of cycles
    n_cycles = len(cycles)
    indices = [i for i in range(n_cycles)]
    # add nodes to the graph
    graph.add_nodes_from(indices)
    # edges of the graph: which cycles are connected
    # by sharing a gauge check of the right type
    cycles_connected: list[tuple[int, int]] = []
    for i, j in itertools.combinations(indices, r=2):
        # get the gauge checks which are in both cycles
        common_gauges = set(cycles[i].gauges) & set(cycles[j].gauges)
        # check if they share a gauge of the cleaned type, if so create an edge
        if any(gauge.type == pauli_to_clean for gauge in common_gauges):
            cycles_connected.append((i, j))
    graph.add_edges_from(cycles_connected)
    # find all cycles that are 'closed' inside the window
    closed_cycles = [k for k, cycle in enumerate(cycles) if cycle.is_closed]
    # initialize the list of cycles to open
    open_cycles: set[int] = set()
    # for each connected component of cycles
    for component in networkx.connected_components(graph):
        # if all cycles are closed then add it to the set of
        # closed components and move on
        if all([k in closed_cycles for k in component]):
            continue
        open_cycles |= component
    return open_cycles


def update_superstabilizers(
    superstabilizers: set[SuperStabilizer], gauges_along_closed_holes: set[Stabilizer]
) -> set[Stabilizer]:
    """Helper function that updates the superstabilizers of the hole by keeping
    only the gauge checks that were found to be around closed cycles in the hole.
    The gauge checks are the open cycles would no longer form or be part of
    superstabilizers.

    Input arguments:
        `superstabilizers`: A set of all the superstabilizers in the open hole.
        `gauges_along_closed_holes`: All the gauge checks found to be around
        closed cycles within the open hole.

    Output arguments:
        A set of superstabilizers.
    """
    # find all remaining superstabilizers and update gauge checks to keep
    # along the boundary (they become stabilizers)
    updated_superstabilizers: set[Stabilizer] = set()
    for superstabilizer in superstabilizers:
        # get all the gauges in the superstabilizer that are
        # part of closed holes
        gauges = set(superstabilizer.gauges) & gauges_along_closed_holes
        # no gauges, move on
        if len(gauges) == 0:
            continue
        # define a graph whose connected components will be the remaining
        # superstabilizers
        gauge_graph = networkx.Graph()
        # for each pair of gauges in the initial superstabilizer
        # (in the extended window), we determine if they should still be
        # combined in the initial window
        for g1, g2 in itertools.combinations(gauges, r=2):
            # to determine if they should be combined, we use all the
            # the gauge checks along closed holes
            for g in gauges_along_closed_holes:
                # if the gauge check is of the same type as the
                # superstabilizer just keep going
                if g.type == superstabilizer.type:
                    continue
                # if both gauges g1 and g2 share data qubits with
                # the bridging gauge g, then add an edge between
                # g1 and g2 in the graph
                g1g_qubits = set(g1.data_qubits) & set(g.data_qubits)
                g2g_qubits = set(g2.data_qubits) & set(g.data_qubits)
                if any(g1g_qubits) and any(g2g_qubits):
                    gauge_graph.add_edge(g1, g2)
                    break
        # each connected component of the graph is a superstabilizer
        # that we add the set of all superstabilizers in the hole
        for component in networkx.connected_components(gauge_graph):
            updated_superstabilizers.add(SuperStabilizer(list(component)))
    return updated_superstabilizers


def clean_gauges_on_edge(
    cycles: list[Cycle],
    window: Window,
    pauli_to_keep: PauliT,
    superstabilizers: set[SuperStabilizer],
) -> BoundaryDeformation:
    """Function that cleans gauges of the wrong type in an open hole along the
    edge of the patch.

    The function works as follows:
    1) We get the position of the hole in the patch, i.e. along which edge
    (left, bottom, right or top) it is.
    2) We find all the cycles in the hole that are 'open', i.e. truncated at
    the boundary.
    3) For each open cycle, we collect which gauges to keep and which to
    gauges to remove along that cycle. We also collect the data qubits
    along that cycle since they will be part of the deformed boundary.
    4) For each closed cycle, we also collect the gauges along that cycle.
    They will form smaller superstabilizers.
    5) We create and return a `BoundaryDeformation` object.

    Input arguments:
        `cycles`: A list of `Cycle` objects associated with all the cycles
        (or contours) within the open hole to clean. (There can be more than one.)
        `window`: A `Window` object containing information about the initial window.
        `pauli_to_keep`: Pauli type of the gauge checks to keep along the boundary.
        `superstabilizers`: A set of all the superstabilizers in the open hole.

    Output arguments:
        A `BoundaryDeformation` object.
    """
    # define the Pauli type to clean along the boundary
    pauli_to_clean = PauliT.X if pauli_to_keep == PauliT.Z else PauliT.Z
    # get the position of the boundary deformation
    pos = get_boundary_position(cycles, window, pauli_to_clean)

    # find open cycles to clean
    open_cycles = find_open_cycles(cycles, pauli_to_clean)
    # instantiate gauges to keep along the boundary
    gauges_to_keep: set[Stabilizer] = set()
    # instantiate gauges to clean along the boundary
    gauges_to_clean: set[Stabilizer] = set()
    # we collect all the data qubits around the open hole
    bdy_str: set[Pos] = set()
    # we keep all gauges around closed components (that will form superstabilizers)
    gauges_along_closed_holes: set[Stabilizer] = set()
    # for each connected component of cycles
    for k, cycle in enumerate(cycles):
        if k in open_cycles:
            # get all the gauge checks that need to be kept or cleaned
            new_gauges = set(gauge for gauge in cycle.gauges)
            gauges_to_keep |= set(gauge for gauge in new_gauges if gauge.type == pauli_to_keep)
            gauges_to_clean |= set(gauge for gauge in new_gauges if gauge.type == pauli_to_clean)
            # update the boundary with data qubits from the open cycles
            bdy_str |= set(q for q in cycle.qubits if window.in_window(q))
            bdy_str |= set(
                q for gauge in new_gauges if gauge.type == pauli_to_clean for q in gauge.data_qubits
            )
        else:
            gauges_along_closed_holes |= set(cycle.gauges)
    gauges_along_closed_holes -= gauges_to_clean

    # update the superstabilizers of the hole now that some gauges are cleaned
    updated_superstabilizers = update_superstabilizers(superstabilizers, gauges_along_closed_holes)
    # remove the gauges of superstabilizers from gauges we keep along boundary
    # (they are already accounted for within the superstabilizers)
    gauges_to_keep -= set(SuperStabilizer.decompose(list(updated_superstabilizers)))

    # we append the boundary deformations
    return BoundaryDeformation(
        updated_superstabilizers,
        gauges_to_keep,
        BoundaryDeformationStrings([BoundaryDeformationString(bdy_str, pos)]),
        [],
    )
