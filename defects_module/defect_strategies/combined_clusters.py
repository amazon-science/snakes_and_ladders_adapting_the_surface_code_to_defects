"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx

from defects_module.base import (
    Pos,
    Stabilizer,
    SuperStabilizer,
)

from .code_distance import compute_effective_distance
from .make_windows import Window


@dataclass
class CombinedCluster:
    """A dataclass that contains all DefectClusters inside a window or patch,
    which are combined together. The class stores the combined superstabilizers,
    the undamaged stabilizers in the patch and the effective distance of the patch.
    This class basically has all the information needed about the strategy being
    applied, and the properties needed to build a patch object.
    """

    window: Window
    """Window of the patch in which live the clusters"""
    effective_distance: tuple[int, int]
    """Effective distance of patch in which the CombinedCluster lives"""
    x_superstabilizers: list[SuperStabilizer]
    """X-type superstabilizers in the window"""
    z_superstabilizers: list[SuperStabilizer]
    """Z-type superstabilizers in the window"""
    undamaged_stabilizers: list[Stabilizer]
    """Undamaged stabilizers in the window"""
    horizontal_distance_graph: networkx.Graph
    """Detector graph used to compute the horizontal effective distance"""
    vertical_distance_graph: networkx.Graph
    """Detector graph used to compute the vertical effective distance"""
    boundaries: dict[str, list[Pos]]
    """Boundaries of the deformed patch"""
    efficient_strategy: list = field(default_factory=list)
    """List of efficient checks in the patch."""
    holes: list = field(default_factory=list)
    """List of Holes in the patch."""

    @property
    def stabilizers(self) -> list[Stabilizer]:
        """Property that returns the stabilizers of the patch."""
        return self.x_superstabilizers + self.z_superstabilizers + self.undamaged_stabilizers

    @property
    def num_data_qubits(self) -> int:
        """Returns the number of active data qubits in the patch."""
        data_qubits = []
        for stabilizer in self.stabilizers:
            data_qubits += stabilizer.data_qubits
        return len(set(data_qubits))

    def __equal__(self, cluster: CombinedCluster) -> bool:
        """Method that checks if this object is equal to another.

        Input arguments:
            `cluster`: A `CombinedCluster` object with which the
            current object will be compared.

        Output arguments:
            True or False.
        """
        return set(self.stabilizers) == set(cluster.stabilizers)

    def __in_list__(self, cluster_list: list[CombinedCluster]) -> bool:
        """Method that checks if this object is equal to any other
        object in a list.

        Input arguments:
            `cluster_list`: A list of `CombinedCluster` objects: each
             one of them will be compared with the current object.

        Output arguments:
            True or False.
        """
        return any([self.__equal__(cluster) for cluster in cluster_list])

    def set_effective_distance(self, window: Window) -> None:
        """Method that computes the effective distance given
        the window in which lives the patch.

        Input arguments:
            `window`: Window object containing information
            about the window in which lives the patch.
        """
        # compute the effective distance
        (
            self.effective_distance,
            self.horizontal_distance_graph,
            self.vertical_distance_graph,
        ) = compute_effective_distance(
            window,
            self.efficient_strategy,
            self.x_superstabilizers,
            self.z_superstabilizers,
            self.boundaries,
        )
