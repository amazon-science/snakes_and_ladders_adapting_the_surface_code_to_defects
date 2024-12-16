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
)
from defects_module.defect_strategies.boundary_deformation.utility import (
    BoundaryDeformation,
    BoundaryDeformationString,
    BoundaryDeformationStrings,
)
from defects_module.defect_strategies.make_windows import Window

from .make_cycles import Cycle
from .utility import InvalidHoleFoundError, any_on_boundary


def find_corners_in_hole(
    window: Window, qubits_on_contours: set[Pos]
) -> list[tuple[tuple[int, int], Pos]]:
    """Helper function that finds which corner falls within the open hole.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `qubits_on_contours`: All the data qubits that are along the contours or
        cycles within the open hole.

    Output arguments:
        A tuple containing
        1) a tuple of integers
        2) a Pos
        Both of them identify the corner in the hole.
    """
    # check if there are any qubits along each boundary
    x0 = any([node.x == window.xlims[0] + 1 for node in qubits_on_contours])
    x1 = any([node.x == window.xlims[1] - 1 for node in qubits_on_contours])
    y0 = any([node.y == window.ylims[0] + 1 for node in qubits_on_contours])
    y1 = any([node.y == window.ylims[1] - 1 for node in qubits_on_contours])
    # instantiate the list of corners in the hole
    corners: list[tuple[tuple[int, int], Pos]] = []
    # check which corner is in the hole
    if y0 and x0:
        # bottom left
        corners.append(((1, 0), Pos(window.xlims[0] + 1, window.ylims[0] + 1)))
    if y1 and x0:
        # top left
        corners.append(((3, 0), Pos(window.xlims[0] + 1, window.ylims[1] - 1)))
    if y0 and x1:
        # bottom right
        corners.append(((1, 2), Pos(window.xlims[1] - 1, window.ylims[0] + 1)))
    if y1 and x1:
        # top right
        corners.append(((3, 2), Pos(window.xlims[1] - 1, window.ylims[1] - 1)))
    # NOTE: we can only do up to two corners per hole at the moment
    if len(corners) > 2:
        raise InvalidHoleFoundError(
            "Strategy dropped:"
            " Found more than two corners in the hole, which is not"
            " currently supported."
        )
    # return the corner
    return corners


def clean_gauges_at_corner(
    cycles: list[Cycle], window: Window, gauges: list[Stabilizer]
) -> list[BoundaryDeformation]:
    """Function that cleans gauge checks of the wrong type in an open corner hole
    for each possible corner position.

    The function works as a follows:
    1) For each cycle within the hole, we check if the cycle touches the boundary.
    If not, we move to the next cycle.
    2) We then check that the new corner is within the initial window. We then remove
    it from the `gauge_graph` of the `Cycle` object such that we can find two connected
    components in that graph (because we removed the corner connecting them). These
    two components correspond to two boundary deformations along two perpendicular edges
    of the patch. If we found more than two connected components, we move to the next cycle.
    3) For these two connected components, we evaluate which one could be of vertical pauli
    type and which could be of horizontal pauli type. They could be interchangeable, depending
    on the hole.
    4) For each gauge checks in `gauges`, we determine if the gauge check falls in the vertical
    or horizontal boundary. We then determine if the gauge needs to be cleaned. If it is cleaned,
    then we also update the data qubits in the new boundary by adding those contained in the
    cleaned gauge check.
    5) For each possibility (each corner, each combination of vertical/horizontal placement) we
    return a `BoundaryDeformation` object.

    Input arguments:
        `cycles`: A list of `Cycle` objects associated with all the cycles
        (or contours) within the open hole to clean. (There can be more than one.)
        `window`: A `Window` object containing information about the initial window.
        `gauges`: List of all the gauge checks within the open hole.

    Output aeguments:
        A list of BoundaryDeformation objects.
    """
    # define the corner type by finding the vertical horizontal position of the corner
    qubits_on_contours = {q for cycle in cycles for q in cycle.qubits}
    initial_corners = find_corners_in_hole(window, qubits_on_contours)

    # The code only supports 1 and 2 corners
    # NOTE: There is no optimization ATM for the two corner case.
    if len(initial_corners) == 2:
        # Case with two corners
        pos1, initial_corner1 = initial_corners[0]
        pos2, initial_corner2 = initial_corners[1]
        if {pos1, pos2} == {(1, 0), (3, 0)} or {pos1, pos2} == {(1, 2), (3, 2)}:
            # vertical
            pauli_to_keep = PauliT.Z if window.vertical_logical == PauliT.X else PauliT.X
            pos = 0 if {pos1, pos2} == {(1, 0), (3, 0)} else 2
        else:
            # horizontal
            pauli_to_keep = PauliT.X if window.vertical_logical == PauliT.X else PauliT.Z
            pos = 1 if {pos1, pos2} == {(1, 0), (1, 2)} else 3
        # determine all gauge checks that need to be kept along the boundary
        # instantiate gauges to keep along the boundary
        gauges_to_keep: set[Stabilizer] = set()
        # for each gauge in the punctured superstabilizers, we determine if we keep
        # it or not and update the data qubits on the deformed boundary string
        gauges = [gauge for cycle in cycles for gauge in cycle.gauges]
        for gauge in gauges:
            if gauge.type == pauli_to_keep:
                # keep the gauge check along the deformed boundary if of the right type
                gauges_to_keep.add(gauge)
        data_qubits = {q for gauge in gauges_to_keep for q in gauge.data_qubits}
        new_corner1, new_corner2 = sorted(
            list(itertools.combinations(data_qubits, r=2)),
            key=lambda x: (x[0] - x[1]).x ** 2 + (x[0] - x[1]).y ** 2,
        )[-1]

        if pauli_to_keep != window.vertical_logical:

            if (initial_corner1.y < initial_corner2.y and new_corner1.y < new_corner2.y) or (
                initial_corner1.y > initial_corner2.y and new_corner1.y > new_corner2.y
            ):
                corners = [
                    (initial_corner1, new_corner1),
                    (initial_corner2, new_corner2),
                ]
            else:
                corners = [
                    (initial_corner1, new_corner2),
                    (initial_corner2, new_corner1),
                ]

        if pauli_to_keep == window.vertical_logical:

            if (initial_corner1.x < initial_corner2.x and new_corner1.x < new_corner2.x) or (
                initial_corner1.x > initial_corner2.x and new_corner1.x > new_corner2.x
            ):
                corners = [
                    (initial_corner1, new_corner1),
                    (initial_corner2, new_corner2),
                ]
            else:
                corners = [
                    (initial_corner1, new_corner2),
                    (initial_corner2, new_corner1),
                ]

        return [
            BoundaryDeformation(
                set(),
                gauges_to_keep,
                BoundaryDeformationStrings(
                    [
                        BoundaryDeformationString(data_qubits, pos),
                    ]
                ),
                corners,
            )
        ]

    # Case with one corner
    (horizontal_pos, vertical_pos), initial_corner = initial_corners[0]
    boundary_deformations: list[BoundaryDeformation] = []

    # test all potential corner positions
    for cycle in cycles:
        if not (
            any_on_boundary(window.xlims, "x", cycle.qubits)
            or any_on_boundary(window.ylims, "y", cycle.qubits)
        ):
            continue
        for corner in cycle.qubits:
            # get the graph of the hole contour
            graph = cycle.gauge_graph.copy()
            if corner not in graph.nodes:
                continue
            # remove the corner
            graph.remove_node(corner)
            # Get the two components corresponding to the vertical and horizontal boundaries
            components = list(networkx.connected_components(graph))
            if len(components) != 2:
                continue
            # Assign which component is on the horizontal and vertical boundaries:
            # there might be more than one possibility! It depends how the hole modifies
            # the patch and which corner we are considering
            horizontal_possibilites = [
                component
                for component in components
                if any_on_boundary(window.ylims, "y", component)
            ]
            vertical_possibilities = [
                component
                for component in components
                if any_on_boundary(window.xlims, "x", component)
            ]
            # This should not be needed, but in case we use a cycle along the edge
            # that is contained in a corner hole, it might happen that we do not
            # find a string along one of the boundaries using the previous approach.
            # In this case, we just pass both strings as possibilities.
            if len(horizontal_possibilites) == 0:
                horizontal_possibilites = [component for component in components]
            if len(vertical_possibilities) == 0:
                vertical_possibilities = [component for component in components]
            # make sure the horizontal and vertical components are unique
            for combination in itertools.product(horizontal_possibilites, vertical_possibilities):
                horizontal, vertical = combination
                if horizontal == vertical:
                    continue
                horizontal.add(corner)
                vertical.add(corner)
                # determine all gauge checks that need to be kept along the boundary
                # instantiate gauges to keep along the boundary
                gauges_to_keep = set()
                # for each gauge in the punctured superstabilizers, we determine if we keep
                # it or not and update the data qubits on the deformed boundary string
                for gauge in gauges:
                    # check if the gauge check is along the vertical boundary string:
                    # it needs to have at least two data qubits along the string
                    # but could have more if the string goes around it
                    in_vertical = len([q for q in gauge.data_qubits if q in vertical]) >= 2
                    # or the horizontal boundary string
                    in_horizontal = len([q for q in gauge.data_qubits if q in horizontal]) >= 2
                    # check if the gauge Pauli type is the same as the vertical
                    # logical Pauli type
                    equal_pauli_types = gauge.type == window.vertical_logical
                    if (in_vertical and not equal_pauli_types) or (
                        in_horizontal and equal_pauli_types
                    ):
                        # keep the gauge check along the deformed boundary if of the right type
                        gauges_to_keep.add(gauge)
                    elif in_vertical and equal_pauli_types:
                        # if not of the right type then update the vertical
                        # string if gauge on vertical
                        vertical |= set(q for q in gauge.data_qubits)
                    elif in_horizontal and not equal_pauli_types:
                        # or update horizontal string if on horizontal
                        horizontal |= set(q for q in gauge.data_qubits)

                # we append the boundary deformations
                boundary_deformation = BoundaryDeformation(
                    set(),
                    gauges_to_keep,
                    BoundaryDeformationStrings(
                        [
                            BoundaryDeformationString(vertical, vertical_pos),
                            BoundaryDeformationString(horizontal, horizontal_pos),
                        ]
                    ),
                    [(initial_corner, corner)],
                )
                boundary_deformations.append(boundary_deformation)
    return boundary_deformations
