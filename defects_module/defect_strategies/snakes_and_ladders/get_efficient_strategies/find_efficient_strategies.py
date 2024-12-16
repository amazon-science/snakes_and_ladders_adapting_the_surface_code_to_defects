"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import copy
import itertools
from typing import Sequence

import networkx

from defects_module.base import (
    Pos,
)

from .apply_efficient_checks import (
    efficient_repurposing_per_cluster,
)
from .define_heuristics import (
    Heuristics,
    Scores,
    decide_if_we_keep_the_strategy,
    filter_clusters,
)
from .efficient_defect_clusters import (
    EfficientDefectsCluster,
)


def _find_efficient_strategies_per_link_assignment(
    cluster: EfficientDefectsCluster,
    scores: Scores,
    repurpose_ancillas: bool = True,
) -> list[EfficientDefectsCluster]:
    """Function that finds maximal ways of reviving zombie data qubits, where
    maximal is defined as using the efficient strategy over the Auger method
    wherever it is possible. Generates all possible combinations of efficient
    checks that are valid, i.e. such that any ancilla is only repurposed once.
    If a heuristic is given, it narrows the search to a subset of the combinations.
    If repurpose_ancillas is set to False, just return the result of Auger's
    algorithm.

    The function works as follows:
        1) If repurpose_ancillas = False, then we run the Auger method. To do so,
        we call the EfficientDefectsCluster.zombify_data_qubit_from_missing_checks
        method to zombify all data qubits around the defective or zombie ancillas.
        We also zombify data qubits in weight-1 checks.
        2) If repurpose_ancillas = True, then we find all possible combinations
        of efficient checks with efficient_repurposing_per_cluster. For each
        possible combination of efficient checks we zombify the data qubits that
        around ancilla defects/zombies that are not part of efficient checks. We
        also zombify all data qubits that are part of weight-1 checks. We decide
        if we keep the strategy or not based on the heuristics.

    Input arguments:
        `cluster`: An EfficientDefectsCluster object that has information about
        the defect configuration.
        `scores`: A Scores objects containing information about the heuristics
        and the optimal configurations.
        `repurpose_ancillas`: A boolean that determines if we use the efficient
        strategy (True) or not (False) in which case we use Auger only. The
        default value is True.

    Output arguments:
        A list of EfficientDefectsCluster objects corresponding to the different
        strategies for the defects cluster.
    """
    if not repurpose_ancillas:
        ### AUGER STRATEGY ###
        # Create a copy of the EfficientDefectsCluster
        updated_cluster = copy.deepcopy(cluster)
        # Evaluate ancilla pair combination: keep only possible efficient checks,
        # enable or disable zombie data qubits and check the validity of the strategy.
        updated_cluster.zombify_data_qubit_from_missing_checks(None)
        # Determine which of the sentenced qubits are zombies and assign them.
        # Here we assume that the relived zombies are now active data qubits.
        updated_cluster.zombify_frozen_qubits()
        return [updated_cluster]

    ### SNAKES N LADDERS STRATEGY ###

    # Get possible pairs of `efficient` checks for each ancilla defect in the cluster
    efficient_checks = efficient_repurposing_per_cluster(cluster)
    # We identify all orientations that can be flipped
    flippable_repurposings = [j for j, checks in enumerate(efficient_checks) if len(checks) > 1]
    # We also define an iterator corresponding to the number of orientations to flip at a
    # given time as will become clear below
    n_to_flip = 0
    # List of EfficientDefectsCluster after applying the efficient strategy on them.
    efficient_clusters: list[EfficientDefectsCluster] = []
    # List of connectivity graphs for the efficient checks for the cluster
    graphs: list[networkx.Graph] = []
    # Loop through all combinations of efficient checks (one for each ancilla defect/zombie)
    while True:
        # We flipped all of them already, we stop.
        if n_to_flip > len(flippable_repurposings):
            break
        # Check if we reached our limit with the heuristics, if so we stop
        if scores.check_if_break():
            break
        # Check if we reached the maximum number of solutions to keep
        if scores.n_buffer_per_cluster and len(efficient_clusters) > scores.n_buffer_per_cluster:
            break
        # We get configurations where n_to_flip orientations were flipped from
        # the initial configuration
        for inds in itertools.combinations(flippable_repurposings, r=n_to_flip):
            # We define the combination of weight-2 checks
            checks_comb = [
                checks[1] if j in inds else checks[0] for j, checks in enumerate(efficient_checks)
            ]
            # Create a copy of the EfficientDefectsCluster: it will be updated with
            # the strategy below
            updated_cluster = copy.deepcopy(cluster)
            # Evaluate the combination: keep only possible efficient checks,
            # enable or disable zombie data qubits and check the validity of the strategy.
            updated_cluster.zombify_data_qubit_from_missing_checks(checks_comb)
            # if the strategy is not valid, move on to the next
            if updated_cluster.efficient_checks_graph is None:
                continue
            # Only consider the strategy if it was not considered before
            if any(
                networkx.utils.graphs_equal(updated_cluster.efficient_checks_graph, graph)
                for graph in graphs
            ):
                continue
            # Add the connectivity graph to the list
            graphs.append(updated_cluster.efficient_checks_graph)
            # Determine which of the sentenced qubits are zombies and assign them.
            # Here we assume that the relived zombies are now active data qubits.
            updated_cluster.zombify_frozen_qubits()
            # Decide if we keep the strategy or not based on heuristics (if any)
            scores, efficient_clusters = decide_if_we_keep_the_strategy(
                scores, updated_cluster, efficient_clusters
            )
            # We check if we reached our limit with heuristics, if so we stop
            if scores.check_if_break():
                break
            # We filter the solutions by the distance loss
            efficient_clusters = filter_clusters(efficient_clusters, scores)
            # we check if have enough solutions already, if so we stop
            if (
                scores.n_buffer_per_cluster
                and len(efficient_clusters) > scores.n_buffer_per_cluster
            ):
                break
        # Update the iterator for the number of orientations to flip
        n_to_flip += 1

    # Return all the solutions for this link assignment
    return efficient_clusters


def find_efficient_strategies_per_cluster(
    cluster: EfficientDefectsCluster,
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
) -> list[EfficientDefectsCluster]:
    """Function that iterates over all subsets of the link defects in the given cluster,
    i.e. for each link defect we disable either the neighboring ancilla qubit or the data
    qubit. For each of these configurations, we apply the Snakes & Ladders algorithm to
    find all valid strategies or reduce the search space according to the given heuristic.
    We return the strategies found for each configuration.

    The function works as follows:
        1) Create all possible assignments for the link defects to data defects or
        zombie ancillas.
        2) For each assignment, we call `_find_efficient_strategies_per_link_assignment`
        above to find all valid strategies.
        3) We filter the strategies to keep the most optimal ones if heuristics were
        passed.

    Input arguments:
        `cluster`: An EfficientDefectsCluster object that has information about
        the defect configuration.
        `repurpose_ancillas`: A boolean that determines if we use the efficient
        strategy (True) or not (False) in which case we use Auger only. The
        default value is True.
        `use_heuristics`: A Heuristics object that contains information about which
        heuristics to use.

    Output arguments:
        A list of EfficientDefectsCluster objects corresponding to the different
        strategies for the defects cluster.
    """
    # find all combinations of data defects and ancilla zombies from link defects
    # make a new EfficientDefectsClusters for each combination
    efficient_clusters: list[EfficientDefectsCluster] = []
    # here we consider only link defects where the data and ancilla are NOT already
    # defective: otherwise, we are already taking care of them by taking care of the
    # data and ancilla defects
    only_link_defects: list[Sequence[Pos]] = [
        link for link in cluster.links_defective if set(link).isdisjoint(cluster.defects)
    ]
    if not repurpose_ancillas:
        # we map link defects to data defects, with Auger's algorithm
        only_link_defects = [
            (node,) for link in only_link_defects for node in link if node in cluster.data_qubits
        ]
    # if the heuristic says so, we only consider mapping a link defect to an ancilla zombie
    elif use_heuristics.link_defect_to_ancilla:
        only_link_defects = [
            (node,) for link in only_link_defects for node in link if node in cluster.ancillas
        ]
    scores = Scores(use_heuristics)
    # for each combination of data defects and ancilla zombies converted from the link defects
    for defects in itertools.product(*only_link_defects):
        data_defects = set(defects) & cluster.data_qubits
        ancilla_zombie = set(defects) & cluster.ancillas_active
        # create a new EfficientDefectsCluster
        new_cluster = copy.deepcopy(cluster)
        # update its fields
        new_cluster.defects |= data_defects
        new_cluster.data_qubits_defective |= data_defects
        new_cluster.data_qubits_zombie -= new_cluster.data_qubits_defective
        new_cluster.data_qubits_active -= data_defects
        new_cluster.ancillas_zombie |= ancilla_zombie
        # Find strategies for this version of the cluster
        efficient_clusters.extend(
            _find_efficient_strategies_per_link_assignment(
                new_cluster,
                scores=scores,
                repurpose_ancillas=repurpose_ancillas,
            )
        )
        scores.partial_reset()

    # Filter the list of efficient defect clusters with heuristics
    efficient_clusters = filter_clusters(efficient_clusters, scores)

    return efficient_clusters
