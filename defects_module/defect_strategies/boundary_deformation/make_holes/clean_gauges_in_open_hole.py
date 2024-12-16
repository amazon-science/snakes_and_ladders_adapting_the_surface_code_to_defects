"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import (
    PauliT,
    SuperStabilizer,
)
from defects_module.defect_strategies.boundary_deformation.utility import (
    BoundaryDeformation,
)
from defects_module.defect_strategies.make_windows import Window

from .clean_gauges_at_corner import (
    clean_gauges_at_corner,
)
from .clean_gauges_on_edge import (
    clean_gauges_on_edge,
)
from .make_cycles import Cycle
from .utility import any_on_boundary


def clean_gauges_if_open(
    cycles: list[Cycle], window: Window, superstabilizers: set[SuperStabilizer]
) -> list[BoundaryDeformation]:
    """Function that returns the gauges in an open hole to be kept along the deformed
    boundary that goes around that hole.

    The function works as follows:
    1) We determine if the hole is at a corner or if it's along an edge of the patch.
    2) If it is a corner hole, we call `clean_gauges_at_corner`.
    3) If it is a corner hole, we also check if the hole touches to the edges of the
    extended window. If it does, then we also call `clean_gauges_on_edge` for both
    Pauli types: this is because it makes sense to define the new corner on both
    edges of the patch.  If it does not, then we call `clean_gauges_on_edge` only
    for the Pauli type for which the new corner makes sense. Otherwise, we would
    end up with more than four corners.
    4) If it is not a corner hole but instead an edge hole, we determine the type
    of the edge on which is the hole. We then call `clean_gauges_on_edge`. Unlike
    the corner hole case, there is always a single BoundaryDeformation object for
    this case.
    5) We return the list of possible BoundaryDeformation objects.

    Input arguments:
        `cycles`: A list of `Cycle` objects associated with all the cycles
        (or contours) within the open hole to clean. (There can be more than one.)
        `window`: A `Window` object containing information about the initial window.
        `superstabilizers`: A set of all the superstabilizers in the open hole.

    Output arguments:
        A list of BoundaryDeformation objects.
    """
    boundary_deformations: list[BoundaryDeformation] = []

    # check if the contour touches the vertical and horizontal boundaries
    qubits_on_contours = {q for cycle in cycles for q in cycle.qubits}
    on_vertical_bdy = any_on_boundary(window.xlims, "x", qubits_on_contours)
    on_horizontal_bdy = any_on_boundary(window.ylims, "y", qubits_on_contours)

    # check if the hole is at a corner
    if on_vertical_bdy and on_horizontal_bdy:

        # Do the case with the most digging (the entire hole is potentially opened)
        gauges = [gauge for superstabilizer in superstabilizers for gauge in superstabilizer.gauges]
        boundary_deformations.extend(clean_gauges_at_corner(cycles, window, gauges))

        # Do the cases where we minimize the amount of digging by putting the
        # corner on an endpoint of smallest opened connected component and treating
        # the corner hole as an edge hole
        for pauli_to_keep in [PauliT.X, PauliT.Z]:
            boundary_deformations.append(
                clean_gauges_on_edge(cycles, window, pauli_to_keep, superstabilizers)
            )

    else:
        # otherwise it's just an edge hole

        # we keep the z-type checks only if we are either
        # (1) on vertical boundary of Z type (X logical on the vertical)
        # (2) on horizontal boundary of Z type (Z logical on the vertical)
        # otherwise we keep the x-type checks
        horizontal_logical = PauliT.Z if window.vertical_logical == PauliT.X else PauliT.X
        pauli_to_keep = window.vertical_logical if on_horizontal_bdy else horizontal_logical

        # Clean the gauges of the wrong type of the open connected component
        boundary_deformations.append(
            clean_gauges_on_edge(cycles, window, pauli_to_keep, superstabilizers)
        )

    return boundary_deformations
