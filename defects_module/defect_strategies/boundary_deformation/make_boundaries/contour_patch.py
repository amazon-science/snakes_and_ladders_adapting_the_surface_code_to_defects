"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import (
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.defect_strategies.boundary_deformation.utility import BoundaryDeformation
from defects_module.defect_strategies.make_windows import Window

from .draw_boundaries import make_boundaries
from .remove_extrusions import remove_extrusions


def contour_patch(
    window: Window,
    boundary_deformation: BoundaryDeformation,
    undamaged_stabilizers: set[Stabilizer],
    superstabilizers: set[SuperStabilizer],
) -> tuple[set[Stabilizer], set[SuperStabilizer], dict[str, list[Pos]]]:
    """Function that returns the boundaries of the patch given the
    boundary defomration and the stabilizers. The stabilizers could be updated
    when trying to draw the boundaries: some stabilizers that fall outside
    the boundaries are removed.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `boundary_deformation`: A `BoundaryDeformation` object containing information
        about all the boundary deformations around the open holes in the patch.
        `undamaged_stabilizers`:  A set of all the undamaged stabilizers in the patch,
        defined in the initial window.
        `superstabilizers`: A set of all the superstabilizers in the patch, also defined
        in the initial window.

    Output arguments:
        A tuple containing:
        1) The updated set of undamaged stabilizers in the patch.
        2) The updated set of superstabilizers in the patch.
        3) A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as lists of Pos objects, containg information about the boundary along
        each of the four sides of the patch. The list of Pos is the list of data qubits used
        to define the boundaries.
    """
    # make the boundaries around the patch
    boundaries = make_boundaries(
        window,
        boundary_deformation,
        undamaged_stabilizers | superstabilizers,
    )
    # Remove extrusions from the patch (all the stabilizers outside the boundaries)
    undamaged_stabilizers, superstabilizers = remove_extrusions(
        boundaries, undamaged_stabilizers, superstabilizers
    )
    return undamaged_stabilizers, superstabilizers, boundaries
