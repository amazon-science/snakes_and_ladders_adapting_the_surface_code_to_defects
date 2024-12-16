"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .efficient_defect_clusters import (
    EfficientDefectsCluster,
)


@dataclass
class Heuristics:
    """Class that defines what heuristics to use to speed up the
    computation of the best strategies. The default is to have no
    heuristics at all.

    Input arguments:
        `link_defects_to_ancilla`: Whether or not we map all link
        defects to ancilla zombies. The default is False. If True,
        then link defects are never mapped to data defects. If False,
        then link defects are mapped to both data defects and ancilla
        zombies.
        `n_zombie`: The number of configurations with the same number
        of functioning data qubits being recorded before we stop the
        search. If a strategy is found to have more functioning data
        qubits then the counter is restarted. This heuristic uses the
        fact that we order the strategies such that efficient checks
        that are not lossy (that locally minimizes the distance loss)
        are tried first. This is a strong heuristic which makes the
        runtime polynomial but can lead to a reduction of the effective
        distance. The default value is None, which means it is not
        applied and keep all valid strategies. n_zombie must be set
        to an integer bigger than 0 if used.
        `n_max`: The number of smallest max distance loss values
        being kept. Keeping more than n_max = 1 values helps with
        the global maximization of the effective distance. The default
        value is None, which means it is not applied and we keep all
        valid strategies. n_max must be set to an integer bigger than
        0 if used.
        `n_sum`: The number of smallest sum distance loss values
        being kept. Keeping more than n_sum = 1 values helps with
        the global maximization of the effective distance. The default
         value is None, which means it is not applied and we keep all
        valid strategies. n_sum must be set to an integer bigger than
        0 if used.
        `n_skip`: The minimum number of strategies to evaluate for each
        assignment of link defects to data defects or ancilla zombies
        before moving to the next assignment. We only move on if we
        tried n_skip samples that were sub-optimal and were not recorded.
        The default value is None, which means it is not applied and we
        keep all valid strategies. n_skip must be set to an integer bigger
        than 0 if used.
    """

    link_defect_to_ancilla: bool = False
    """A strong heuristic that restricts the search to the
    strategies that only map link defects to ancilla defects."""
    n_zombie: Optional[int] = None
    """Number of configurations to keep that have the same
    number of zombie qubits, before stopping the search for
    each assignment of link defects.
    """
    n_max: Optional[int] = None
    """Number of minimal max distance loss values to keep. Default
    value is None which corresponds to keeping all."""
    n_sum: Optional[int] = None
    """Number of minimal sum distance loss values to keep. Default
    value is None which corresponds to keeping all."""
    n_skip: Optional[int] = None
    """
    Number of sub-optimal strategies to evaluate for each assignment
    of link defects before moving on to the next assignement. Here
    sub-optimal means that the strategy was not worth recording because
    the distance drop was greater than the ones previously recorded.
    """
    n_sol_max_per_cluster: Optional[int] = None
    """
    Maximum number of strategies to keep per cluster. If set to 1 then
    there is no exponential build-up in memory.
    """
    n_buffer_per_cluster: Optional[int] = None
    """
    Maximum number of strategies to keep per cluster during the search.
    The default is `n_sol_max_per_cluster`. This allows the user to
    increase the search space while keeping only a few final strategies.
    """
    n_sol_max: Optional[int] = None
    """Maximum number of strategies allowed for the patch.
    """

    def __post_init__(self):
        if self.n_zombie is not None:
            assert isinstance(self.n_zombie, int) and self.n_zombie > 0
        if self.n_max is not None:
            assert isinstance(self.n_max, int) and self.n_max > 0
        if self.n_sum is not None:
            assert isinstance(self.n_sum, int) and self.n_sum > 0
        if self.n_skip is not None:
            assert isinstance(self.n_skip, int) and self.n_skip > 0
        if self.n_sol_max_per_cluster is not None:
            assert isinstance(self.n_sol_max_per_cluster, int) and self.n_sol_max_per_cluster > 0
            if self.n_buffer_per_cluster is not None:
                assert (
                    isinstance(self.n_buffer_per_cluster, int)
                    and self.n_buffer_per_cluster >= self.n_sol_max_per_cluster
                )
            else:
                self.n_buffer_per_cluster = self.n_sol_max_per_cluster
        if self.n_sol_max is not None:
            assert isinstance(self.n_sol_max, int) and self.n_sol_max > 0


@dataclass
class StrategyCounter:
    """Class that helps count the number of strategies for
    a given cluster.

    Input arguments:
        `n_skip`: The minimum number of strategies to evaluate
        for each assignment of link defects to data defects or
        ancilla zombies before moving to the next assignment.
    """

    n_skip: Optional[int] = None
    """Number of sub-optimal strategies before switching to
    the next link to data defect/ancilla zombie assignment.
    """
    n_eval: int = 0
    """Number of strategies evaluated so far in this
    assignment.
    """
    n_recorded: int = 0
    """Number of strategies that were recorded."""
    break_flag: bool = False
    """Flag that determined if the search keeps going."""

    def record(self, recorded: bool) -> None:
        """Method that records if a strategy was kept or not.

        It updates the number of evaluated strategies n_eval
        by 1, and the number of recorded strategies n_recorded
        by 1 if `recorded` is True. It also checks if we need
        to stop the search: if the number of evaluated strategies
        n_eval is bigger or equal to the allowed number n_skip,
        then we set break_flag to True indicating that we move
        to the next link defect assignment.

        Input arguments:
            `recorded`: True if the strategy is kept, False otherwise.
        """
        if self.n_skip:
            self.n_eval += 1
            if recorded:
                self.n_recorded += 1
            if self.n_recorded == 0 and self.n_eval >= self.n_skip:
                self.break_flag = True

    def reset(self) -> None:
        """Method that resets the internal counters."""
        self.n_eval = 0
        self.n_recorded = 0
        self.break_flag = False


@dataclass
class ZombieCounter:
    """Class that helps filter the strategies for a given cluster
    using the number of data qubits that are alive, i.e. by
    minimizing the number of zombie data qubits.

    Input arguments:
        `n_zombie`: The number of configurations with the same number
        of functioning data qubits being recorded before we stop the
        search.
    """

    n_zombie: Optional[int] = None
    """Maximum number of strategies with the same number of
    functioning data qubits that is maximal.
    """
    n_active: int = 0
    """Highest number of functioning data qubits recorded."""
    n_configs: int = 0
    """Number of strategies with the same number of functioning
    data qubits.
    """
    break_flag: bool = False
    """Flag that determines if the search keeps going or not."""

    def record(self, n_active: int) -> bool:
        """Method that updates the internal counters.

        It resets the counter n_configs to 0 if the number
        of active qubits n_active is larger than the
        maximum number that was recorded so far. Otherwise,
        if it is equal than it updates n_configs by 1.
        We then check if n_configs is greater than the allowed
        number n_zombie: if so, we set break_flag to True,
        indicating that we stop the search. For both of the
        cases where n_active is bigger or equal, we keep the
        strategy and return True. If n_active is instead smaller
        than we update nothing and return False.

        Input arguments:
            `n_active`: The number of functioning data
            qubits in the strategy.

        Output arguments:
            True if we keep the strategy, else False.
        """
        if self.n_zombie is None:
            # we do not filter and thus record the new strategy
            # (we just don't track it internally)
            return True
        # we cannot add more strategies
        if self.break_flag:
            return False
        # the strategy has fewer functional data qubits than the best one before
        if n_active < self.n_active:
            return False
        # if this strategy has the same number of functional data qubits as the
        # best one we evaluated before
        if n_active == self.n_active:
            self.n_configs += 1
        # check if we reached the maximum number of configurations
        if self.n_configs >= self.n_zombie:
            self.break_flag = True
        return True

    def update(self, n_active: int) -> None:
        if self.n_zombie is None:
            return
        # if this strategy has more functional data qubits (aka fewer zombies) than
        # all strategies evaluated before
        if n_active > self.n_active:
            self.n_configs = 0
            self.n_active = n_active

    def reset(self) -> None:
        """Method that resets the internal counters."""
        self.n_active = 0
        self.n_configs = 0
        self.break_flag = False


@dataclass
class DistanceLossFilter:
    """Class that helps filter the strategies for a given cluster using
    the distance loss (max or sum).

    Input arguments:
        `filter_fcn`: A string ('max' or 'sum') that defines the function
        to apply on the distance loss of the defects cluster (which is a
        tuple of integers) and returns a single integer. The default value
        is None, meaning that this heuristic/filtering is not applied.
        `n_val`: The number of minimial distance loss values to keep.
        The default value is None, when this heuristic is not applied.
        If used, then n_val must be an integer bigger than 0.
    """

    filter_fcn: Optional[str] = None
    """Minimize the max distance loss (between horizontal and vertical)
    per cluster."""
    n_val: Optional[int] = None
    """Keep a strategy only if its max distance loss is among the
    `number_of_max_distance_loss_values` lowest values (not counting
    duplicates)."""
    # set the list of values
    values: list[int] = field(default_factory=list)
    """Store the values of max distance loss that we have already encountered"""

    def __post_init__(self) -> None:
        # if we don't set the number of values to keep then
        # we also set the filter function to None (no heuristic)
        if self.n_val is None:
            self.filter_fcn = None
        # make sure the filter function is max or sum if not None
        # and that we have an integer number of strategies > 0
        if self.filter_fcn:
            assert isinstance(self.n_val, int) and self.n_val > 0
            assert self.filter_fcn in ["max", "sum"]

    def record(self, distance_loss: tuple[int, int]) -> tuple[bool, list[int]]:
        """Method that updates the best distance loss values that
        were recorded so far.

        It checks if the new distance loss is better than the ones
        previously recorded. If not, nothing is updated and we
        return False, meaning the strategy is not kept. If it is
        better then we add it to the internal list `values` if the
        number of values is less than `n_val`, otherwise we update
        the last element of the list. We then sort the list. We
        return True in this case because the strategy is kept.

        Input arguments:
            `distance_loss`: A tuple of integers corresponding
            to the code distance loss.

        Output arguments:
            The returned boolean indicates whether the new
            strategy will be kept and the list of new values
            if the strategy is to be recorded.
        """
        if self.filter_fcn is None:
            # we do not filter so we return True to record
            # the new strategy (we just don't track it internally)
            return True, []
        assert isinstance(self.n_val, int)
        # compute the distance drop using max or sum depending on the heuristic
        new_value = max(distance_loss) if self.filter_fcn == "max" else sum(distance_loss)
        values = self.values.copy()
        if len(values) == self.n_val and new_value > values[-1]:
            # drop the strategy if we already have n_val different distance
            # drop values and the new strategy is worse than all
            return False, values
        if len(values) < self.n_val:
            # if we haven't recorded n_val different distance drop values,
            # add the new value
            values.append(new_value)
        elif new_value not in values:
            # if we already hit the maximum number of values,
            # and the new value is unique, we get rid of
            # the worse distance currently in the list
            values[-1] = new_value
        # we sort all the values and make sure they are unique
        values = sorted([d for d in set(values)])
        # we return True since we updated the values
        # and will record the new strategy
        return True, values

    def update(self, values: list[int]) -> None:
        """Update the internal values."""
        if self.filter_fcn is None:
            return
        self.values = values

    def filter_clusters(
        self, clusters: list[EfficientDefectsCluster]
    ) -> list[EfficientDefectsCluster]:
        """Method that filters a list of efficient clusters based on
        their distance loss.

        Input arguments:
            `clusters`: A list of EfficientDefectsCluster objects
            corresponding to all strategies that were recorded so far.

        Output arguments:
            A list of all strategies that are kept because their
            distance loss is less or equal to the worst distance loss
            that we wanted to keep.
        """
        if self.filter_fcn is None:
            # we do not do any filtering
            return clusters

        def fcn(x):
            return max(x) if self.filter_fcn == "max" else sum(x)

        # keep a strategy if its max/sum distance loss is among the n_val lowest values
        filtered_clusters = [c for c in clusters if fcn(c.distance_loss) <= self.values[-1]]
        return filtered_clusters


class Scores:
    """Class that helps apply the heuristic to strategies per cluster and keeps track
    of optimal strategies that were found before.

    The constructor works as follows:
        1) We define a StrategyCounter object, used for the heuristic that
        skips to the next link defect assignment to data defect / zombie ancilla
        if no improvement is noted in the distance loss.
        2) We define a ZombieCounter object, used for the heuristic that
        counts the number of zombie data qubits in the cluster and minimizes it.
        3) We define a DistanceLossFilter object for the max distance loss, used
        to filter the strategies by the max of their distance loss in both directions
        such that we keep a fixed number of strategies per cluster.
        4) Same as 3) but for the sum distance loss instead.

    Input arguments:
        `heuristics`: A Heuristics object that contains all the information
        necessary to define the heuristics.
    """

    def __init__(self, heuristics: Heuristics) -> None:
        # heuristic to minimize the number of strategies we go through
        self.skip_to_next_assignment = StrategyCounter(n_skip=heuristics.n_skip)
        # heuristic to minimize the number of zombie data qubits
        self.minimize_zombie_qubits = ZombieCounter(n_zombie=heuristics.n_zombie)
        self.compute_zombie_qubits = bool(heuristics.n_zombie)
        # heuristics to minimize the distance loss
        self.minimize_max_distance_loss = DistanceLossFilter(
            filter_fcn="max", n_val=heuristics.n_max
        )
        self.minimize_sum_distance_loss = DistanceLossFilter(
            filter_fcn="sum", n_val=heuristics.n_sum
        )
        # heuristic to keep a max number of strategies per cluster
        self.n_sol_max_per_cluster = heuristics.n_sol_max_per_cluster
        self.n_buffer_per_cluster = heuristics.n_buffer_per_cluster

    def check_if_break(self) -> bool:
        """Method that adjusts the break flag if we need to skip to the next link
        defect assignment because we either didn't find:
            1) any better strategies with the `skip_to_next_assignment` heuristic
            2) any strategies with more active data qubits with the
            `minimize_zombie_qubits` heuristic.
        """
        # if we have already found enough configurations with the min number of zombies
        # we will stop looking at more configurations
        return self.skip_to_next_assignment.break_flag or self.minimize_zombie_qubits.break_flag

    def partial_reset(self) -> None:
        """Method that partially resets the Scores object before applying to the
        strategies from another link defect assignment. Reset the counters for
        `minimize_zombie_qubits` to allow more exploration. Also reset the counters
        for `skip_to_the_next_assignment`.
        """
        self.skip_to_next_assignment.reset()
        self.minimize_zombie_qubits.reset()


def decide_if_we_keep_the_strategy(
    scores: Scores,
    updated_cluster: EfficientDefectsCluster,
    updated_clusters: list[EfficientDefectsCluster],
) -> tuple[Scores, list[EfficientDefectsCluster]]:
    """Function that is called for each combination of efficient checks (for the same assignment
    of link defects to data defects / ancilla zombies). It keeps track of the strategies that have
    already been checked, and tells the user if the current strategy should be kept or not.

    The function works as follows:
        1) If the `minimize_zombie_qubits` heuristic is used, then we update its internal
        counters and determine if the strategy is worse than the previous ones. If so,
        we don't keep the strategy, otherwise we keep going.
        2) If the `minimize_max_distance_loss` heuristic is used, then we update its
        internal values and determine if the strategy is worth keeping. If it is worth
        keeping then we keep going.
        3) We do the same as 2) but for `minimize_sum_distance_loss`.
        4) We update the `skip_to_next_assignement` heuristic and update the list of strategies
        if we decided to keep the current strategy.

    Input arguments:
        `scores`: a Scores object that records the information needed for applying heuristics.
        `updated_cluster`: The current strategy being considered.
        `updated_clusters`: The list of strategies that have already been considered and recorded.

    Output arguments:
        A tuple containing:
        1) The updated Scores object.
        2) The updated list of strategies.
    """
    # decide if we record the strategy or not
    record: bool = True
    # We compute the distance loss and check for invalid
    # configurations
    distance_loss = updated_cluster.compute_distance_loss()
    updated_cluster.distance_loss = distance_loss
    # if we found invalid superstabilizers when computing
    # the distance loss, we do not keep the strategy.
    if min(distance_loss) == -1:
        return scores, updated_clusters
    # If we use heuristics based on the number of zombie (data) qubits
    if scores.compute_zombie_qubits:
        # number of functional data qubits in the cluster
        n_active = len(updated_cluster.data_qubits_active)
        record = scores.minimize_zombie_qubits.record(n_active)
    # If we use heuristics based on the distance drop per cluster
    if record:
        # Minimize max distance drop
        record, max_values = scores.minimize_max_distance_loss.record(distance_loss)
    if record:
        # Minimize sum distance drop
        record, sum_values = scores.minimize_sum_distance_loss.record(distance_loss)
    # Append the updated cluster if we decided to keep the strategy
    if record:
        if scores.compute_zombie_qubits:
            scores.minimize_zombie_qubits.update(n_active)
        scores.minimize_max_distance_loss.update(max_values)
        scores.minimize_sum_distance_loss.update(sum_values)
        updated_clusters.append(updated_cluster)
    # Decide if we should skip to the next link assignment
    scores.skip_to_next_assignment.record(record)
    return scores, updated_clusters


def filter_clusters(
    updated_clusters: list[EfficientDefectsCluster],
    scores: Scores,
) -> list[EfficientDefectsCluster]:
    """Function that filters a list of defect clusters (the same cluster with different
    assignment of link defects to ancilla and data defects, and/or with different strategies
    applied) based on distance loss in the cluster. This removes the inferior strategies that
    are recorded before we evaluate better ones.

    Input arguments:
        `updated_clusters`: A list of EfficientDefectsCluster objects, i.e. all strategies
        that were kept so far.
        `scores`: A Scores object that contain information about each heuristic and what
        are currently the best values for the distance loss.

    Output arguments:
        A list of EfficientDefectsCluster objects corresponding to the best strategies
        so far in terms of distance loss.
    """
    if len(updated_clusters) <= 1:
        return updated_clusters
    # filter the strategies based on the maximum distance loss across two directions
    # keep a strategy if its max distance loss is among the n_max lowest values
    updated_clusters = scores.minimize_max_distance_loss.filter_clusters(updated_clusters)
    # filter the strategies based on the sum of the distance losses across two directions
    # keep a strategy if its total distance loss is among the n_sum lowest values
    updated_clusters = scores.minimize_sum_distance_loss.filter_clusters(updated_clusters)
    # restrain the number of strategies per cluster
    if isinstance(scores.n_sol_max_per_cluster, int) and len(updated_clusters) > 0:
        # we sort the unique distance loss values
        distance_losses = list({c.distance_loss for c in updated_clusters})

        def sorting_func(distance_loss):
            return max(distance_loss), sum(distance_loss)

        distance_losses.sort(key=sorting_func)
        # we partition the max number strategies per unique distance loss value
        counts = {j: 0 for j, _ in enumerate(distance_losses)}
        n = scores.n_buffer_per_cluster
        assert isinstance(n, int)
        m = len(distance_losses)
        for k in range(m):
            for j in range(m - k):
                counts[j] += n // (m - k)
            n -= (n // (m - k)) * (m - k)
        # we keep the best unique distance loss strategies (while minimizing zombie qubits)
        _updated_clusters: list[EfficientDefectsCluster] = []
        for j, count in counts.items():
            _clusters = [c for c in updated_clusters if c.distance_loss == distance_losses[j]]
            _clusters.sort(key=lambda x: x.data_qubits_active, reverse=True)
            _updated_clusters += _clusters[0:count]
            if j < len(counts) - 1:
                counts[j + 1] += len(_clusters) - count
        updated_clusters = _updated_clusters
    # return the filtered strategies
    return updated_clusters
