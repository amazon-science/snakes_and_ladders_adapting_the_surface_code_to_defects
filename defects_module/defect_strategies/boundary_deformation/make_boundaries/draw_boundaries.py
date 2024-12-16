"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools

import networkx

from defects_module.base import (
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.defect_strategies.boundary_deformation.utility import BoundaryDeformation
from defects_module.defect_strategies.make_windows import Window


def _initial_boundaries(window: Window) -> list[set[Pos]]:
    """Function that returns the boundaries of the undamaged patch as a list
    of sets of qubits for each boundary (left, bottom, right, top).

    Input arguments:
        `window`: A `Window` object containing information about the initial window.

    Output arguments:
        A list of sets of qubits in Pos format.
    """
    # We define the boundaries with the data qubits (no ordering)
    # 0: left, 1: bottom, 2: right, 3: top.
    boundaries: list[set[Pos]] = []
    # first we add all the data qubits on the initial boundaries that
    # are still part of the patch
    for k in range(4):
        if k % 2 == 0:
            # left or right
            bdy = set(
                Pos(window.xlims[k // 2] + 1 - k, y)
                for y in range(window.ylims[0] + 1, window.ylims[1] + 1, 2)
            )
        else:
            # bottom or top
            bdy = set(
                Pos(x, window.ylims[(k - 1) // 2] + 2 - k)
                for x in range(window.xlims[0] + 1, window.xlims[1] + 1, 2)
            )
        boundaries.append(bdy)
    return boundaries


def _make_graphs(
    window: Window, boundary_deformation: BoundaryDeformation, stabilizers: set[Stabilizer]
) -> dict[int, networkx.Graph]:
    """Function that creates the graph for each boundary (left, bottom, right, top)
    used to define the boundary (as a a path in the graph). The graphs are returned as
    a dictionary whose keys are the positions of the boundaries, i.e.
    (0: left, 1: bottom, 2: right, 3: top) and whose values are the networkx graphs.

    Each graph is built from the connectivity of the stabilizers and gauges in the
    superstabilizers: each edge connects two data qubits if they are connected to the
    same ancilla qubit. The graph nodes are only data qubits that are either
    part of the initial non-deformed boundaries OR part of the boundary deformations
    obtained from the open holes.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `boundary_deformation`: A `BoundaryDeformation` object containing information
        about all the boundary deformations around the open holes in the patch.
        `stabilizers`: A set of all the stabilizers in the patch, defined in the
        initial window.

    Output arguments:
        A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as networkx graphs.
    """
    # We define the boundaries with the data qubits (no ordering)
    boundaries = _initial_boundaries(window)
    # we define the data qubits in the patch
    data_qubits = set(qubit for stab in stabilizers for qubit in stab.data_qubits)
    boundaries = [bdy & data_qubits for bdy in boundaries]
    # then we add all the data qubits that are along the open holes
    for string in boundary_deformation.boundary_strings.strings:
        boundaries[string.position] |= string.data_qubits & data_qubits
    # we define the graph for each boundary
    graphs = {k: networkx.Graph() for k in range(4)}
    for stab in SuperStabilizer.decompose(list(stabilizers)):
        for k in [0, 2]:
            # determine in which graph the edges should be added
            ind = k + int(stab.type == window.vertical_logical)
            # add the edges
            qubits_in_stab = [q for q in stab.data_qubits if q in boundaries[ind]]
            graphs[ind].add_edges_from(list(itertools.combinations(qubits_in_stab, r=2)))
    # return graphs
    return graphs


class NoBoundariesFoundError(Exception):
    """Error that is raised when a data qubit or an ancilla qubit is on the boundary
    of a patch or a window. The code is then considered to be uncorrectable.
    """

    pass


def _boundaries_from_graphs(graphs: dict[int, networkx.Graph]) -> list[list[Pos]]:
    """Functions that identifies boundaries from the gauge graph by first identifying
    possible corners, and second find the shortest path between each pair of those
    corners. We keep the shortest paths overall.

    Input arguments:
        `graphs`: A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as networkx graphs.

    Ouput arguments:
        A list of list of Pos, corresponding to data qubits along the left, right,
        bottom and top boundaries.
    """
    # update the boundaries with the qubits that are in the stabilizers
    # (might be different because we cleaned some frozen qubits)
    boundaries: list[set[Pos]] = [set(graphs[k].nodes) for k in range(4)]
    # we define the potential corner positions
    corners_prod = itertools.product(
        boundaries[3] & boundaries[0],
        boundaries[0] & boundaries[1],
        boundaries[1] & boundaries[2],
        boundaries[2] & boundaries[3],
    )
    corners_list = [
        list(corners) + [corners[0]] for corners in corners_prod if len(set(corners)) == 4
    ]
    if not corners_list:
        # failure to draw boundaries
        raise NoBoundariesFoundError(
            "Strategy dropped: Did not find four corners to draw boundaries."
        )
    # reorder the qubits in the boundaries
    # we find the shortest path between the corners
    ordered_boundaries_list: list[list[list[Pos]]] = []
    for corners in corners_list:
        valid_corners = True
        paths: list[list[Pos]] = []
        # we find a path for each side of the patch
        for k in range(4):
            # we find the source and target nodes
            # from the chosen set of corners
            source, target = corners[k : k + 2]
            try:
                # we attempt to find the shortest path
                paths.append(networkx.shortest_path(graphs[k], source=source, target=target))
            except Exception:
                # if it fails, we deem the corners to be invalid
                valid_corners = False
                break
        # if we found valid corners, we append the paths (or boundaries) to the list
        if valid_corners:
            ordered_boundaries_list.append(paths)
    # If no boundaries were found, we raise an error.
    if not ordered_boundaries_list:
        raise NoBoundariesFoundError(
            "Strategy dropped: Did not find paths for all four boundaries."
        )
    # We take the boundaries that have the shortest length overall (closer to rectangle)
    ordered_boundaries = sorted(
        ordered_boundaries_list, key=lambda x: len([i for j in x for i in j])
    )[0]
    # we return the boundaries
    return ordered_boundaries


def _count_checks(window: Window, qubit: Pos, stabilizers: set[Stabilizer]) -> list[int]:
    """Helper function that computes the number of X and Z checks applied on each data qubit.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `qubit`: A Pos corresponding to the data qubit being inspected.
        `stabilizers`: A set of all the stabilizers in the patch.

    Output arguments:
        A list containing the number of X and Z checks.
    """
    # we define a list to count how many X and Z checks are being
    # applied on the data qubit
    checks = [0, 0]
    # we iterative over all the stabilizers in the patch
    for stab in stabilizers:
        # if the stabilizer is not acting on it then we skip it
        if qubit not in stab.data_qubits:
            continue
        # we always order the list as number of checks of
        # horizontal, vertical pauli type
        checks[stab.type == window.vertical_logical] += 1
    return checks


def _validate_boundaries(
    window: Window, stabilizers: set[Stabilizer], ordered_boundaries: list[list[Pos]]
) -> None:
    """Function that validates a choice of boundaries. This function only raises an error
    if the boundaries are invalid.

    We check the following:
    1) No opposite boundaries must be touching.
    2) Every pair of adjescent boundaries share a single data qubit.
    3) There are four corners.
    4) The number of X and Z checks on each data qubit is correct.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `stabilizers`: A set of all the stabilizers in the patch.
        `ordered_boundaries`: A list of list of Pos, corresponding to data qubits
        along the left, right, bottom and top boundaries.
    """
    # first we make sure that the logicals meet at all four corners
    # otherwise it is not a valid configuration
    # first we check that opposite boundaries do not cross
    # otherwise the configuration is invalid
    if any(set(ordered_boundaries[0]) & set(ordered_boundaries[2])) or any(
        set(ordered_boundaries[1]) & set(ordered_boundaries[3])
    ):
        raise NoBoundariesFoundError("Strategy dropped: Had opposite boundaries touching.")
    # second we count the number of corners and make sure there are exactly 4
    final_corners: set[Pos] = set()
    for bdy1, bdy2 in itertools.combinations(ordered_boundaries, r=2):
        new_corner = set(bdy1) & set(bdy2)
        # if ever the two boundaries cross more than once
        # then it's not a valid configuration
        if len(new_corner) > 1:
            raise NoBoundariesFoundError(
                "Strategy dropped: Had more than one qubit shared between two boundaries."
            )
        final_corners |= new_corner
    if len(final_corners) != 4:
        raise NoBoundariesFoundError("Strategy dropped: Found less than four corners.")
    # third we get all the data qubits that are along the boundary
    qubits_on_bdy = set(
        ordered_boundaries[0]
        + ordered_boundaries[1]
        + ordered_boundaries[2]
        + ordered_boundaries[3]
    )
    # and we iterate over the stabilizers in the patch
    # for each data qubit on the boundary we check which stabilizers
    # are acting on it
    for qubit in qubits_on_bdy:
        # count the number of checks on the data qubit
        checks_count = _count_checks(window, qubit, stabilizers)
        # get the target number of checks
        if qubit in final_corners:
            # case where the qubit is a corner
            target_count = [1, 1]
        elif qubit in ordered_boundaries[0] + ordered_boundaries[2]:
            # case where the qubit is along the vertical
            # (therefore we need to have 2 checks of the horizontal pauli type)
            target_count = [2, 1]
        else:
            # case where the qubit is along the horizontal
            # (therefore we need to have 2 checks of the vertical pauli type)
            target_count = [1, 2]
        # we check if the number of X and Z checks on that data qubit is valid
        if checks_count != target_count:
            raise NoBoundariesFoundError(
                f"Strategy dropped: Had frozen qubit {qubit} along boundary "
                "(probably due to an invalid corner placement)."
            )


def make_boundaries(
    window: Window, boundary_deformation: BoundaryDeformation, stabilizers: set[Stabilizer]
) -> dict[str, list[Pos]]:
    """Function that creates the boundaries of the patch by returning a dictionary
    with the list of qubits for each boundary (left, bottom, right, top).

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `boundary_deformation`: A `BoundaryDeformation` object containing information
        about all the boundary deformations around the open holes in the patch.
        `stabilizers`: A set of all the stabilizers in the patch.
        `ordered_boundaries`: A list of list of Pos, corresponding to data qubits
        along the left, right, bottom and top boundaries.

    Output arguments:
        A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as lists of Pos objects, containg information about the boundary along
        each of the four sides of the patch. The list of Pos is the list of data qubits used
        to define the boundaries.
    """
    # we make graphs from the stabilizers to define the boundaries
    graphs = _make_graphs(window, boundary_deformation, stabilizers)
    # extract the boundaries from the graphs
    ordered_boundaries = _boundaries_from_graphs(graphs)
    # make sure the boundaries are valid
    _validate_boundaries(window, stabilizers, ordered_boundaries)
    # Now that the boundaries are valid, we build the dictionary to be returned
    labels = ["left", "bottom", "right", "top"]
    return dict([item for item in zip(labels, ordered_boundaries, strict=True)])
