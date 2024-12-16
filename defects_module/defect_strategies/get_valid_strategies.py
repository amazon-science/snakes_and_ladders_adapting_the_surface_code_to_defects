"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools
import warnings
from typing import Sequence

from defects_module.base import (
    Pos,
)
from defects_module.code_library import (
    RotatedSurfaceCode,
)

from .boundary_deformation import run_boundary_deformation
from .combined_clusters import CombinedCluster
from .defects_clusters import generate_defects_clusters
from .make_windows import Window
from .snakes_and_ladders import Heuristics, run_snakes_and_ladders


class NSolMaxReachedError(Exception):
    """Error that is raised when too many strategies where found
    in the extended window before the boundary deformation. Note
    that the total runtime is not necessarily strongly correlated
    with this number but the memory is. In general, if more
    strategies are found the higher computer resources are needed.
    The total runtime also depends on the different possibilities
    for the boundary deformation, particularly for corner holes.
    """

    pass


def get_valid_strategies_in_patch(
    patch: RotatedSurfaceCode,
    ancilla_defects: set[Pos],
    data_defects: set[Pos],
    link_defects: set[tuple[Pos, Pos]],
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
    add_padding: bool = True,
) -> list[CombinedCluster]:
    """Function that returns all the possible strategies for the defective surface code patch.
    These strategies include different choices of superstabilizers and different possible boundary
    deformations, which were all checked to be valid.

    The function works as follows:
    1. We define the *extended* and *initial* windows.
    2. We make sure that the defects configuration is valid. It happens at the second step
    because the initial window includes padding ancillas.
    3. We separate code ancilla defects from padding ancilla defects.
    4. We generate the defect clusters in the *extended* window.
    5. We run the Snakes and Ladders algorithm (if repurpose_ancillas = True, else Auger)
    in the *extended* window using the defect clusters found in 4. At this step, we use
    heuristics if any were passed. By default, runtime and memory grow exponentially
    in the number of ancilla and link defects (~ 2^n). Different heuristics can make
    this growth sub-exponential or polynomial instead.
    6. For each possible set of superstabilizers (strategy) we run the boundary deformation
    algorithm in the *initial* window: for each hole found in the *extended* window, we need
    to determine if it is closed or open in the *initial* window. For each possible boundary
    deformation, we collect all the information necessary to build a patch (i.e. a
    `CombinedCluster` object) and compute the effective code distance.
    7. We return a list of all possible strategies, i.e. a list of `CombinedCluster` objects.
    Each CombinedCluster object has fields for the effective distance, superstabilizers and
    undamaged stabilizers.

    Input arguments:
        `patch`: A RotatedSurfaceCode object which is the rotated surface code
        patch in which are found the defective components.
        `ancilla_defects`: A set of all the ancilla defects in the patch.
        `data_defects`: A set of all the data defects in the patch.
        `link_defects`: A set of all the link defects in the patch.
        `repurpose_ancillas`: An option to use snakes and ladders instead of Auger.
        The default is True.
        `use_heuristics`: Option to add heuristics when finding the different possible
        superstabilizers with snakes and ladders.
        `add_padding`: An option to add padding ancillas.

    Output arguments:
        A list of `CombinedCluster` objects corresponding each to a unique, valid strategy.
    """
    # make the initial and extended windows
    window, extended_window = Window.make_windows(patch, add_padding=add_padding)

    # make sure that the defect configuration is valid
    all_window_ancillas = window.ancillas | window.padding_ancillas
    assert ancilla_defects.issubset(all_window_ancillas)
    assert data_defects.issubset(window.data_qubits)
    for link in link_defects:
        assert isinstance(link, Sequence) and len(link) == 2
        assert link[0] - link[1] in [Pos(1, 1), Pos(1, -1), Pos(-1, 1), Pos(-1, -1)]
        assert (link[0] in window.data_qubits and link[1] in all_window_ancillas) or (
            link[1] in window.data_qubits and link[0] in all_window_ancillas
        )

    # separate the code defects from the padding defects:
    # We keep track of the defective qubits in the padding only to know where it is impossible to
    # repurpose the ancillas. These defects are otherwise not perceived as 'defects'
    # in snakes & ladders, i.e. they are only used when we run efficient_repurposing
    unavailable_padding_ancillas = ancilla_defects & window.padding_ancillas
    updated_ancilla_defects = ancilla_defects - unavailable_padding_ancillas
    updated_link_defects = link_defects.copy()
    # for each defective link that has a padding qubit on one end, we mark the padding qubit
    # as unavailable but don't consider the link as a defect that needs to be handled
    for link in link_defects:
        if set(link).isdisjoint(window.padding_ancillas):
            continue
        updated_link_defects.remove(link)
        unavailable_padding_ancillas |= set(link) & window.padding_ancillas

    # find all DefectClusters (i.e. clusters of defects) in the window.
    clusters = generate_defects_clusters(
        extended_window,
        updated_ancilla_defects,
        data_defects,
        updated_link_defects,
        unavailable_padding_ancillas,
    )
    print(f"Found {len(clusters)} defect clusters.")

    # List of all DefectClusters after using the efficient strategy
    efficient_clusters: list[list] = []
    # EfficientDefectsCluster as DefectsCluster
    # For each DefectCluster in the window, 'generate_efficient_clusters' finds where
    # we can use the efficient strategy and returns the superstabilizers. Each DefectCluster
    # can be solved with different deformations of the stabilizers depending on where the
    # efficient strategy is used, stored in 'updated_clusters'. We consider all of them in
    # 'efficient_clusters' because the total effective distance of the window is impacted by
    # how the local deformations connect in the error strings.
    for cluster in clusters:
        updated_clusters = run_snakes_and_ladders(
            cluster,
            repurpose_ancillas=repurpose_ancillas,
            use_heuristics=use_heuristics,
        )
        if len(updated_clusters) == 0 and repurpose_ancillas:
            print(
                cluster.ancillas_defective, cluster.data_qubits_defective, cluster.links_defective
            )
            warnings.warn(
                "Snakes and Ladders failed. Trying Auger.",
                UserWarning,
                stacklevel=2,
            )
            updated_clusters = run_snakes_and_ladders(
                cluster,
                repurpose_ancillas=False,
                use_heuristics=use_heuristics,
            )
        if len(updated_clusters) == 0:
            # if one of the defect clusters has no strategy, then we did not suceed in
            # finding an overall strategy
            warnings.warn(
                "Did not find any strategies for one or more clusters.",
                UserWarning,
                stacklevel=2,
            )
            return []
        efficient_clusters.append(updated_clusters)

    # Get an estimate of the total distance for each possible
    # combination of local strategies per cluster
    combinations = list(
        itertools.product(*[list(range(len(strategies))) for strategies in efficient_clusters])
    )
    # Sort all possible combinations with respect to their distance loss
    def get_total_distance_loss(combination):
        distance_loss = [0, 0]
        for i, j in enumerate(combination):
            dist = efficient_clusters[i][j].distance_loss
            # we do not store the distance loss for Auger,
            # so just use dist = (0, 0) for that case
            if not dist:
                continue
            distance_loss[0] += dist[0]
            distance_loss[1] += dist[1]
        return max(distance_loss), sum(distance_loss)

    combinations.sort(key=get_total_distance_loss)
    # Restrict to maximum number of strategies (we take the ones that we think
    # will lead to minimal distance loss across the patch)
    if use_heuristics.n_sol_max:
        combinations = combinations[0 : use_heuristics.n_sol_max]

    # List of all combined clusters in the window
    combined_clusters: list[CombinedCluster] = []
    for combination in combinations:
        cluster_combination = [efficient_clusters[i][j] for i, j in enumerate(combination)]
        # we determine which holes are open and which are closed in
        # the initial window, and fix the gauges around the open holes
        # we then get the stabilizers given the current strategy with fixed
        # closed and open holes in the window: we have different combinations
        # because some open holes have more than one valid strategies (i.e corner placement)
        _combined_clusters = run_boundary_deformation(window, cluster_combination)
        if not _combined_clusters:
            continue

        # for each valid configuration we compute the effective distance
        for _cb in _combined_clusters:
            # make sure it is a new strategy, otherwise move on
            if _cb.__in_list__(combined_clusters):
                continue
            # compute the effective distance
            _cb.set_effective_distance(window)
            # add it to the list of valid strategies for this window
            combined_clusters.append(_cb)

    return combined_clusters
