"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from typing import Sequence

from defects_module.base import (
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.defect_strategies.combined_clusters import CombinedCluster
from defects_module.defect_strategies.defects_clusters import DefectsCluster
from defects_module.defect_strategies.make_windows import Window
from defects_module.defect_strategies.snakes_and_ladders import TriangularCheck

from .get_patch_stabilizers import (
    NoPatchFoundError,
    collect_data_from_holes,
    update_patch_stabilizers,
)
from .make_boundaries import NoBoundariesFoundError, contour_patch
from .make_holes import Hole, InvalidHoleFoundError


def make_strategy(
    window: Window,
    efficient_strategy: list[TriangularCheck],
    holes: list[Hole],
    undamaged_stabilizers: set[Stabilizer],
    superstabilizers: set[SuperStabilizer],
    boundaries: dict[str, list[Pos]],
) -> CombinedCluster:
    """Function that creates a CombinedCluster object in the *initial window*
    given the set of stabilizers in the window and the boundaries. The CombinedCluster
    object will contain all the information necessary to build the patch for this
    strategy.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `undamaged_stabilizers`: A set of all the undamaged stabilizers in the patch,
        defined in the initial window.
        `superstabilizers`: A set of all the superstabilizers in the patch, also define
        in the initial window.
        `boundaries`: A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as lists of Pos objects, containg information about the boundary along
        each of the four sides of the patch. The list of Pos is the list of data qubits used
        to define the boundaries.

    Output arguments:
        A CombinedCluster object.
    """
    # We separate the X and Z type superstabilizers
    x_superstabilizers = [ss for ss in superstabilizers if ss.type == PauliT.X]
    z_superstabilizers = [ss for ss in superstabilizers if ss.type == PauliT.Z]
    # We create a CombinedCluster using the previously computed strategy
    return CombinedCluster(
        window,
        (0, 0),
        sorted(x_superstabilizers, key=lambda s: min(s.ancilla)),
        sorted(z_superstabilizers, key=lambda s: min(s.ancilla)),
        sorted(list(undamaged_stabilizers), key=lambda s: min(s.ancilla)),
        None,
        None,
        boundaries,
        efficient_strategy,
        holes,
    )


def get_strategies_from_holes(
    window: Window, efficient_strategy: list[TriangularCheck], holes: list[Hole]
) -> list[CombinedCluster]:
    """Function that returns a list of `CombinedCluster` objects: one for each possible
    boundary deformation in the initial window given the efficient strategy in the extended
    window. There can be multiple choices for the boundary deformation due to corner placement
    when an open hole includes a corner. Each CombinedCluster object contains all the information
    necessary to build a patch.

    Input arguments:
        `window`: Window object containing information about the *initial window*.
        `holes`: List of `Hole` objects: each object contains information about each
        hole in the patch.

    Output arguments:
        A list of `CombinedCluster` objects, one for each possible strategy
        associated with a specific choice of boundary deformation. The CombinedCluster
        object gives all the information necessary to build the patch.
    """
    stabilizers, boundary_deformations = collect_data_from_holes(window, holes)
    # We define the list of all possible strategies (CombineCluster objects)
    # given the different possible boundary deformations from corner placement
    combined_clusters: list[CombinedCluster] = []
    # We create a strategy for each possible corner placement (if it is valid)
    for boundary_deformation in boundary_deformations:
        # We update the stabilizers given the boundary deformation
        try:
            updated_undamaged_stabilizers, updated_superstabilizers = update_patch_stabilizers(
                stabilizers, boundary_deformation
            )
        except Exception as exception:
            # this task can fail if there are no stabilizer found in the patch
            # after cleaning -- see the function `update_patch_stabilizers` for
            # more information about this exception
            if isinstance(exception, NoPatchFoundError):
                print(exception)
                continue
            else:
                raise exception
        # We make the boundaries around the patch
        try:
            (updated_undamaged_stabilizers, updated_superstabilizers, boundaries,) = contour_patch(
                window,
                boundary_deformation,
                updated_undamaged_stabilizers,
                updated_superstabilizers,
            )
        except Exception as exception:
            # this task can fail for different reasons including an invalid corner placement --
            # see `contour_patch` for more information
            if isinstance(exception, NoBoundariesFoundError):
                print(exception)
                continue
            else:
                raise exception
        # We add the strategy to the list
        combined_clusters.append(
            make_strategy(
                window,
                efficient_strategy,
                holes,
                updated_undamaged_stabilizers,
                updated_superstabilizers,
                boundaries,
            )
        )
    # We return the list of strategies
    return combined_clusters


def run_boundary_deformation(
    window: Window, cluster_combination: Sequence[DefectsCluster]
) -> list[CombinedCluster]:
    """Function that runs the whole algorithm for the boundary deformation.

    The function works as follows:
        1) We determine which holes are open and which are closed in the initial window.
        We then clean the gauges around the open holes.
        2) We then get the stabilizers given the current strategy (i.e. superstabilizers
        defined in the extended window) for each possible boundary deformation: we have
        open holes might have more than one valid deformation due to corner placement.

    Input arguments:
        `window`: Window object containing information about the *initial window*.
        `cluster_combination`: Sequence of `DefectsCluster` objects from which we can
        determine the superstabilizers and the holes. These DefectsCluster objects
        are defined in the *extended* window.

    Output arguments:
        A tuple containing:
        1) A list of `TriangularCheck` objects (efficient checks) in the patch
        2) A list of `CombinedCluster` objects, one for each possible strategy
        associated with a specific choice of boundary deformation. The CombinedCluster
        object gives all the information necessary to build the patch.
    """
    # we determine which holes are open and which are closed in
    # the initial window, and fix the gauges around the open holes
    try:
        (
            efficient_strategy,
            superstabilizers,
            data_defects,
            ancilla_defects,
            ancilla_zombies,
        ) = Hole.collect_data_to_make_holes(cluster_combination)
        holes = Hole.make_holes(
            window, superstabilizers, data_defects, ancilla_defects, ancilla_zombies
        )
    except Exception as exception:
        if isinstance(exception, InvalidHoleFoundError):
            print(exception)
            # sometimes a strategy can be deemed invalid for different reasons
            # like non-commuting stabilizers, overlapping superchecks, etc.
            # and in such a case we move on to the next
            return []
        else:
            raise exception
    # we get the stabilizers given the current strategy with fixed
    # closed and open holes in the window: we have different combinations
    # because some open holes have more than one valid strategies (i.e corner placement)
    combined_clusters = get_strategies_from_holes(window, efficient_strategy, holes)
    # Return the list of efficient checks and a list of all the possible strategies
    # (one for each possible boundary deformation)
    return combined_clusters
