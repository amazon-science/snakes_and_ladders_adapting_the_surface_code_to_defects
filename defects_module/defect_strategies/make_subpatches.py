"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from defects_module.base import SuperStabilizer

from .boundary_deformation import Hole, get_strategies_from_holes
from .combined_clusters import CombinedCluster
from .filter_strategies_by_distance import DistanceFilterFunctions, FilterDistance
from .make_windows import Window
from .snakes_and_ladders import TriangularCheck


def fit_strategy_to_window(
    window: Window, efficient_strategy: list[TriangularCheck], holes: list[Hole]
) -> tuple[list[TriangularCheck], list[Hole]]:
    """Helper function that keeps only the efficient checks and holes that are
    contained inside the new window.

    Input arguments:
        `window`: Window object containing information about the subpatch.
        `efficient_strategy`: A list of efficient (TriangularCheck) checks.
        `holes`: A list of Hole objects.

    Output arguments:
        A tuple containing
        1) A list of efficient (TriangularCheck) checks in the new window
        2) A list of Hole objects in the new window.
    """
    # efficient checks in the subpatch
    new_efficient_strategy = [
        check
        for check in efficient_strategy
        if {check.data1, check.data2}.issubset(window.data_qubits)
    ]

    # when making a subpatch, there is a special case we need to consider:
    # given that there are now actual defects outside the initial window,
    # we can sometimes end up truncating a hole such that only the boundary
    # checks along one boundary were gauge checks of a superstabilizer in the
    # code patch we are trying to make a subpatch in. This means that,
    # we could have a hole whose 'initial' (aka punctured) superstabilizers
    # in this new window that all commute and therefore the hole is considered
    # closed (as it should be) but with these boundary gauges wrongly combined
    # into a superstabilizer. We therefore need to modify the 'initial'
    # superstabilizers of this hole such each boundary check is its own stabilizer.
    # Importantly, we need to make sure to keep only boundary checks that are
    # repurposed.
    repurposed_ancillas = {check.ancilla for check in new_efficient_strategy}

    def check_boundary_superstabilizers(hole: Hole) -> Hole | None:
        # this is only a problem if the hole is closed
        if not hole.is_closed:
            return hole
        stabs = {s for s in hole.superstabilizers["initial"] if len(s.pauli) > 0}
        # this is only a problem if there is only one single Pauli type
        # (meaning the gauge checks all sit on the boundary)
        if len({s.type for s in stabs}) == 2:
            return hole
        gauges = {gauge for s in stabs for gauge in s.gauges}
        gauges = {gauge for gauge in gauges if gauge.only_ancilla in repurposed_ancillas}
        # if there is no repurposed check, this means that there is actually
        # no 'hole' in the subpatch!
        if not gauges:
            return None
        hole.superstabilizers["initial"] = {SuperStabilizer([gauge]) for gauge in gauges}
        return hole

    # holes in the subpatch
    subpatch_qubits = window.data_qubits | window.ancillas | window.padding_ancillas
    new_holes: list[Hole] = []
    for hole in holes:
        defects = hole.data_qubits_defective | hole.data_qubits_zombie
        defects |= hole.ancilla_qubits_defective | hole.ancilla_qubits_zombie
        if not any(defects & subpatch_qubits):
            continue
        new_hole = hole.update(window)
        checked_hole = check_boundary_superstabilizers(new_hole)
        if checked_hole:
            new_holes.append(checked_hole)

    return new_efficient_strategy, new_holes


def find_subpatch_strategies(
    window: Window, efficient_strategy: list[TriangularCheck], holes: list[Hole]
) -> list[CombinedCluster]:
    """Function that returns a list of strategies given a window which
    fits the subpatch, a list of Hole objects in the subpatch and a list
    of efficient checks. The code distance is computed for each strategy.

    Input arguments:
        `window`: Window object containing information about the subpatch.
        `efficient_strategy`: A list of efficient (TriangularCheck) checks
        contained in the subpatch.
        `holes`: A list of Hole objects contained in the subpatch.

    Output arguments:
        A list of CombinedCluster objects, each corresponding to a single
        strategy for the subpatch.
    """
    # keep only the efficient checks and holes in the new window
    efficient_strategy, holes = fit_strategy_to_window(window, efficient_strategy, holes)
    # get the strategies for the new window and sets of efficient checks and holes
    _combined_clusters = get_strategies_from_holes(window, efficient_strategy, holes)
    combined_clusters: list[CombinedCluster] = []
    for _cb in _combined_clusters:
        # make sure it is a new strategy, otherwise move on
        if _cb.__in_list__(combined_clusters):
            continue
        # compute the effective distance
        _cb.set_effective_distance(window)
        # add it to the list of valid strategies for this window
        combined_clusters.append(_cb)
    return combined_clusters


def find_best_subpatch_strategy(
    window: Window, efficient_strategy: list[TriangularCheck], holes: list[Hole]
) -> CombinedCluster | None:
    """Function that returns a the best strategy given a window which
    fits the subpatch, a list of Hole objects in the subpatch and a list
    of efficient checks. The code distance is computed for each strategy.

    Input arguments:
        `window`: Window object containing information about the subpatch.
        `efficient_strategy`: A list of efficient (TriangularCheck) checks
        contained in the subpatch.
        `holes`: A list of Hole objects contained in the subpatch.

    Output arguments:
        A CombinedCluster object corresponding to the best strategy
        strategy for the subpatch in terms of code distance.
    """
    strategies = find_subpatch_strategies(window, efficient_strategy, holes)

    # number of strategies found
    num_strategies = len(strategies)

    if num_strategies == 0:
        return None

    # find the best subpatch
    filter_distance = FilterDistance(filter_func=DistanceFilterFunctions.return_best_strategies)
    strategy = filter_distance.filter(strategies)[0]

    return strategy
