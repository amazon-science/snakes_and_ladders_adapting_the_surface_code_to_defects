"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import SuperStabilizer
from defects_module.defect_strategies.defects_clusters import DefectsCluster

from .get_efficient_strategies import (
    EfficientDefectsCluster,
    Heuristics,
    find_efficient_strategies_per_cluster,
)
from .make_superstabilizers import (
    combine_checks,
    define_checks,
)


def set_superstabilizers(cluster: EfficientDefectsCluster) -> None:
    """Function that find and assigns the superstabilizers of a
    defect cluster.

    The function works as follows:
        1) Define all the gauges in the defects cluster with the
        function `define_checks`.
        2) Combine the gauges into superstabilizers with the
        function `combine_checks`.
        3) Store the superstabilizers in `cluster`.

    Input arguments:
        `cluster`: A EfficientDefectsCluster object which contains
        all the information regarding a specific strategy such as
        the list of efficient checks and the defects in the cluster.
    """
    # Define efficient and gauge checks
    x_checks, z_checks = define_checks(cluster)
    # Combine checks to form superstabilizers
    x_superstabilizers, z_superstabilizers = combine_checks(cluster, x_checks, z_checks)
    # Assign superstabilizers to the cluster
    cluster.x_superstabilizers = x_superstabilizers
    cluster.z_superstabilizers = z_superstabilizers


def _check_valid_superstabilizer(superstabilizer: SuperStabilizer) -> bool:
    """Check if a superstabilizer has a data qubit appeared more than once
    in its gauge checks. If so, return False: it is not considered a valid
    superstabilizer.
    """
    for s1 in superstabilizer.gauges:
        for s2 in superstabilizer.gauges:
            # consider two different gauge checks
            if s1 == s2:
                continue
            # check if they share data qubits
            if any(set(s1.data_qubits) & set(s2.data_qubits)):
                # if so, return False
                return False
    # otherwise the superstabilizer is valid
    return True


def check_valid_superstabilizers(updated_cluster: EfficientDefectsCluster) -> bool:
    """Check if the superstabilizers of an updated defects cluster
    (with a strategy) are all valid.
    """
    for superstabilizer in updated_cluster.x_superstabilizers + updated_cluster.z_superstabilizers:
        if not _check_valid_superstabilizer(superstabilizer):
            return False
    return True


def run_snakes_and_ladders(
    defect_cluster: DefectsCluster,
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
) -> list[EfficientDefectsCluster]:
    """Function that searches and returns a set of valid superstabilizers for a given
    defect cluster. The function will reduce the search space using any heuristic that
    is passed. If repurpose_ancillas is set to True the function will use efficient
    checks (weight-2 checks using repurposed ancillas), otherwise it will run Auger.

    The function works as follows:
        1) We call the function `find_efficient_strategies_per_cluster` from the
        `get_efficient_strategies` module to find all possible strategies.
        Refer to the README in this module for more details on the algorithm.
        2) We then filter the strategies such as to keep only the unique ones. To do
        so, we compare the list of efficient checks of each strategy, and their set of
        superstabilizers.

    Input arguments:
        `defect_cluster`: A defect cluster object defined in the *extended window*.
        There is no notion of boundary deformation in this module. Any defect cluster
        is defined in a window such that holes defined by the Auger strategy will be
        completely closed.
        `repurpose_ancillas`: Boolean which is set to True by default such that we
        allow efficient checks. If set to False, we only use the Auger strategy.
        `use_heuristics`: A Heuristics object which contains the information about
        which heuristics to use for the search of valid strategies.

    Output arguments:
        A list of EfficientDefectsCluster objects where each object contains information
        about a valid strategy (e.g. superstabilizers) for the defects cluster.
    """
    # define an efficient defects cluster object
    cluster = EfficientDefectsCluster(defect_cluster)

    # List of efficient clusters
    # Find all efficient strategies for the cluster, return a list of the same defect
    # cluster updated with the efficient strategy.
    _updated_clusters = find_efficient_strategies_per_cluster(
        cluster, repurpose_ancillas=repurpose_ancillas, use_heuristics=use_heuristics
    )

    # For each updated cluster, find and assign the superstabilizers.
    updated_clusters: list[EfficientDefectsCluster] = []

    for updated_cluster in _updated_clusters:
        # make sure it's a new strategy, otherwise keep going
        # NB: if two clusters have the same triangular checks THEY will have
        # the same distance loss, no matter if data qubits are defective or
        # zombie, or more generally if some other properties of the defects
        # cluster are different
        if any([updated_cluster.__equal_strategies__(c) for c in updated_clusters]):
            continue
        # we set the superstabilizers of the cluster: we also check if they
        # are valid, if not we keep going
        set_superstabilizers(updated_cluster)
        if not check_valid_superstabilizers(updated_cluster):
            continue
        # make sure it is a new strategy, otherwise keep going
        # NB: same as above, two clusters can have the same superstabilizers
        # at this point despite having non-equal efficient checks because of
        # link defects: an asymmetric check can pass as a weight-2 check in
        # a Auger hole
        if any([updated_cluster.__equal_superstabilizers__(c) for c in updated_clusters]):
            continue
        # add the strategy
        updated_clusters.append(updated_cluster)

        # if the user decides to keep only n_sol_max_per_cluster, then stop
        # when we reach the desired number
        if use_heuristics.n_sol_max_per_cluster:
            if len(updated_clusters) == use_heuristics.n_sol_max_per_cluster:
                break

    return updated_clusters
