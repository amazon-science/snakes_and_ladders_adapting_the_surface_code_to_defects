"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from .combined_clusters import (
    CombinedCluster,
)


class DistanceFilterFunctions:
    """Class that holds functions that sort and filter the strategies."""

    @staticmethod
    def basic_ordering(
        strategies: list[CombinedCluster],
    ) -> dict[tuple[int, int], list[CombinedCluster]]:
        """Function that sorts strategies based on (min distance, total distance),
        where the min and the total are both taken over the code distances along the
        two dimensions. It returns a dictionary that groups the strategies by the
        effective code distances.

        Input arguments:
            `strategies`: A list of strategies (`CombinedCluster` objects) to be
            sorted by their effective code distance.

        Output arguments:
            A dictionary where the keys are the tuples (min, sum) effective code
            distance and the values are lists of all strategies (`CombinedCluster`
            objects) that have the same exact effective code distance.
        """
        unique_keys = {
            (min(strategy.effective_distance), sum(strategy.effective_distance))
            for strategy in strategies
        }
        ordered_strategies: dict[tuple[int, int], list[CombinedCluster]] = {
            unique_key: [] for unique_key in unique_keys
        }
        for strategy in strategies:
            ordered_strategies[
                (min(strategy.effective_distance), sum(strategy.effective_distance))
            ].append(strategy)
        for key in ordered_strategies.keys():
            ordered_strategies[key] = sorted(
                ordered_strategies[key], key=lambda x: x.num_data_qubits, reverse=True
            )
        return ordered_strategies

    @staticmethod
    def return_all_strategies(
        ordered_strategies: dict[tuple[int, int], list[CombinedCluster]]
    ) -> list[CombinedCluster]:
        """Function that returns all the strategies in a flattened list, ordered
        by their effective code distance, starting that those having the best min
        distance (and then the best total distance).

        Input arguments:
            `ordered_strategies`: Dictionary returned by `basic_ordering` where the
            keys are the tuples (min, sum) effective code distance and the values are
            lists of all strategies (`CombinedCluster` objects) that have the same
            exact effective code distance.

        Output arguments:
            A list of `CombinedCluster` objects that all share the same effective code
            distance, which was found to be the most optimal one.
        """
        return [
            strategy
            for _, strategies in sorted(ordered_strategies.items(), reverse=True)
            for strategy in strategies
        ]

    @staticmethod
    def return_best_strategies(
        ordered_strategies: dict[tuple[int, int], list[CombinedCluster]]
    ) -> list[CombinedCluster]:
        """Function that returns the strategies in a flattened list, with the same best
        effective code distance following the rule of having the best min distance first,
        followed by the best total distance. There could be more than single strategy.

        Input arguments:
            `ordered_strategies`: Dictionary returned by `basic_ordering` where the
            keys are tuples (min, sum) effective code distance and the values are
            lists of all strategies (`CombinedCluster` objects) that have the same
            exact effective code distance.

        Output arguments:
            A list of `CombinedCluster` objects that all share the same effective code
            distance, which was found to be the most optimal one.
        """
        best_min_distance = max([key[0] for key in ordered_strategies.keys()])
        best_total_distance = max(
            [key[1] for key in ordered_strategies.keys() if key[0] == best_min_distance]
        )
        return [
            strategy
            for key, strategies in ordered_strategies.items()
            if key == (best_min_distance, best_total_distance)
            for strategy in strategies
        ]


class FilterDistance:
    """Class that filters and/or sorts the strategies based on the effective distance.
    User can specify a filter function that uses the output of basic_ordering,
    i.e. that takes (min_distance, total_distance) as input. Built-in
    options are return_best_strategies and return_all_strategies.
    """

    def __init__(self, filter_func=DistanceFilterFunctions.return_best_strategies) -> None:
        """Constructor.

        Input arguments:
            `filter_func`: Function to use for the filtering of the strategies, based
            on their effective code distance. The default value is `return_best_strategies`.
        """
        self.filter_func = filter_func

    def filter(self, strategies: list[CombinedCluster]) -> list[CombinedCluster]:
        """Method that filters the strategies using the `basic_ordering` function
        in the`DistanceFilterFunctions` class and the filter function that was
        passed as input in the constructor.

        Input arguments:
            `strategies`: List of all the strategies to be filtered, defined as
            `CombinedCluster` objects.

        Output arguments:
            A list of `CombinedCluster` objects, sorted by their effective distance.
        """
        ordered_strategies = DistanceFilterFunctions.basic_ordering(strategies)
        return self.filter_func(ordered_strategies)
