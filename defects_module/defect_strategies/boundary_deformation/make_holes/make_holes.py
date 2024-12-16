"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools
from typing import Sequence

import networkx
from shapely.geometry import Point

from defects_module.base import (
    Pauli,
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.defect_strategies.boundary_deformation.utility import (
    BoundaryDeformation,
)
from defects_module.defect_strategies.defects_clusters import DefectsCluster
from defects_module.defect_strategies.make_windows import Window
from defects_module.defect_strategies.snakes_and_ladders import TriangularCheck
from defects_module.defect_strategies.utility import stabilizers_commute

from .clean_gauges_in_open_hole import clean_gauges_if_open
from .make_cycles import Cycle
from .utility import InvalidHoleFoundError


def pair_superstabilizers(
    x_superstabilizers: list[SuperStabilizer], z_superstabilizers: list[SuperStabilizer]
) -> list[tuple[list[SuperStabilizer], list[SuperStabilizer]]]:
    """Helper function that pairs superstabilizers associated with the same hole.

    The function works as follows:
    1) We build a networkx graph whose nodes are the superstabilizers.
    2) For each pair of X and Z superstabilizers we check if their gauge checks
    commute or not. If any pair of gauge checks if found to not commute, then the
    pair of superstabilizers is part of the same hole and is therefored `paired`.
    3) If the pair of superstabilizers if paired, then we add an edge between
    the two superstabilizers in the graph.
    4) We find all the connected components in the graph. Each connected component
    contains all the superstabilizers that are in the same `hole`.
    We return a list of tuples (X superstabilizers, Z superstabilizers) where each
    element in the kth tuple is a list of X or Z superstabilizers contained in the
    kth hole.

    Input arguments:
        `x_superstabilizers`: List of all X superstabilizers in the patch.
        `z_superstabilizers`: List of all Z superstabilizers in the patch.

    Output arguments:
        A list of tuples: each tuple contains two lists. The first is a list of
        X superstabilizers and the second is a list of Z superstabilizers. The kth
        tuple contains the X and Z superstabilizers contained in the kth hole.
    """
    # instantiate a graph whose edges are the paired superstabilizers
    # with non-commuting gauges
    superstabilizer_graph = networkx.Graph()
    # look at all possible pairs of the x and z superstabilizers
    for x_superstabilizer, z_superstabilizer in itertools.product(
        x_superstabilizers, z_superstabilizers
    ):
        # get the list of gauges from the two superstabilizers
        x_gauges = x_superstabilizer.gauges
        z_gauges = z_superstabilizer.gauges
        paired = False
        for x_gauge, z_gauge in itertools.product(x_gauges, z_gauges):
            # check if the two gauges commute
            if not stabilizers_commute(x_gauge, z_gauge):
                # if not then the two superstabilizers are paired
                paired = True
                break
        # if paired, add the edge to the graph
        if paired:
            superstabilizer_graph.add_edge(x_superstabilizer, z_superstabilizer)

    # instantiate a list of paired superstabilizers with non-commuting gauges
    superstabilizer_pairs: list[tuple[list[SuperStabilizer], list[SuperStabilizer]]] = []
    for component in networkx.connected_components(superstabilizer_graph):
        # we combine the superstabilizers whose don't commute with a
        # superstabilizer of the other type
        x_stab = [stab for stab in component if stab.type == PauliT.X]
        z_stab = [stab for stab in component if stab.type == PauliT.Z]
        # and append to the list of pairs
        superstabilizer_pairs.append((x_stab, z_stab))

    return superstabilizer_pairs


class Hole:
    def __init__(
        self,
        window: Window,
        x_superstabilizers: list[SuperStabilizer],
        z_superstabilizers: list[SuperStabilizer],
        data_defects: set[Pos],
        ancilla_defects: set[Pos],
        ancilla_zombies: set[Pos],
    ) -> None:
        """Class that stores information about a given hole in the defective patch.
        Each hole corresponds to a pair of superstabilizers that are defined such as
        to `patch` the defective components. A hole that is closed in the initial window
        (i.e. can be patched with Auger or efficient strategies) will be kept as it is in
        the patch. However, an open hole, due to the defects being on or close to the
        boundary, will deform the patch boundary. The checks that are kept along the
        deformed boundary are computed here.

        The constructor works as follows:
        1) We compute the modified superstabilizers in the *initial window*.
        2) We check if the hole's superstabilizers are still `connected` (i.e touching)
        in the *initial window*.
        3) We check if the modified superstabilizers still commute or not.
        4) The hole is said open if the superstabilizers are no longer `connected` or
        if they don't commute anymore.
        5) We find all the `contours` (closed cycles in a graph whose connectivity is
        determined by the gates in the superstabilizers) in the hole. They help us
        compute the possible boundary deformations.
        6) We determine which defects are contained inside the hole.
        7) If the hole is open then we break the superstabilizers into smaller
        superstabilizers and stabilizers (previously gauge checks) and keep only the
        checks along the `contours` of the hole (only the contours that are found to
        be open in the *inital window*) that match the type of the boundary. We refer
        to this step as `cleaning the boundary`.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
            `x_superstabilizers`: List of X superstabilizers associated with the hole.
            `z_superstabilizers`: List of Z superstabilizers associated with the hole.
            `data_defects`: List of data defects in the patch.
            `ancilla_defects`: List of ancilla defects in the patch.
            `ancilla_zombies`: List of ancilla zombies in the patch.
        """
        # define the superstabilizers of the hole in the extended window
        assert len(x_superstabilizers) != 0 and len(z_superstabilizers) != 0
        self.superstabilizers: dict[str, set[SuperStabilizer]] = {
            "extended": set(x_superstabilizers + z_superstabilizers)
        }
        self.superstabilizers["initial"] = {
            self.puncture_superstabilizer(window, s) for s in self.superstabilizers["extended"]
        }
        """Set of superstabilizers in the hole, in both the extended and initial windows."""

        # check that the superstabilizers are still `connected`, i.e.
        # that their individual gauges don't commute together
        self.is_connected = self.connected(window)
        """Checks if the gauges form a closed contour around the hole
        in the original window.
        """
        # check that the hole is closed by verifying the commutation
        # relation of the punctured superstabilizers which must
        # also be connected
        self.is_closed = self.commute() & self.is_connected
        """Checks if the hole is closed in the original window. """
        # define the contour of the hole in the extended window
        self.set_contours(data_defects, ancilla_defects, ancilla_zombies)
        # filter the qubits in the hole per type
        self.filter_qubits_inside(window, data_defects, ancilla_defects)
        # instantiate the list of boundary deformations
        self.boundary_deformations: list[BoundaryDeformation] = []
        """List of possible boundary deformations around the hole."""
        # `clean` the gauges in the hole if open such that we
        # keep the right gauge checks along the deformed
        # boundary (which goes around the open hole)
        if not self.is_closed:
            self.clean_gauges_if_open(window)

    def update(self, window: Window) -> Hole:
        """Method that "restricts" the hole to the newly given window, thereby possibly
        opening the hole in the process.

        Input arguments:
            `window`: A `Window` object containing information about the new window.

        Output arguments:
            A Hole object.
        """
        x_superstabilizers = [s for s in self.superstabilizers["extended"] if s.type == PauliT.X]
        z_superstabilizers = [s for s in self.superstabilizers["extended"] if s.type == PauliT.Z]
        return Hole(
            window,
            x_superstabilizers,
            z_superstabilizers,
            self.data_qubits_defective,
            self.ancilla_qubits_defective,
            self.ancilla_qubits_zombie,
        )

    def puncture_superstabilizer(
        self,
        window: Window,
        superstabilizer: SuperStabilizer,
    ) -> SuperStabilizer:
        """Method that punctures a superstabilizer by truncating the window from
        the `extended` back to the `initial` size.

        Inpurt arguments:
            `window`: A `Window` object containing information about the initial window.
            `superstabilizer`: The superstabilizer to puncture or modify.

        Output arguments:
            A modified superstabilizer.
        """
        gauges = []
        for gauge in superstabilizer.gauges:
            if gauge.only_ancilla not in window.ancillas | window.padding_ancillas:
                continue
            if gauge.only_ancilla in window.padding_ancillas:
                # When we repurpose an ancilla using two opposite neighboring ancillas,
                # the stabilizers checked by the other pair of opposite neighboring ancillas
                # (perpendicular ancillas) become gauge checks. We need to check that we can
                # form a valid superstabilizer in this way. When the efficient checks are along
                # the boundary such that these perpendicular ancillas include a padding ancilla,
                # we cannot form a superstabilizer because there is no check measured by the padding
                # ancilla in the original code. This special case is handled here.
                # The way we handle this is to make sure that every active padding ancilla is only
                # ever repurposed (i.e. do an efficient check)
                if (
                    gauge.type == window.vertical_logical and gauge.only_ancilla.x in window.xlims
                ) or (
                    gauge.type != window.vertical_logical and gauge.only_ancilla.y in window.ylims
                ):
                    continue
            paulis = [
                Pauli(qubit, PauliT(gauge.type))
                for qubit in gauge.data_qubits
                if window.in_window(qubit)
            ]
            if paulis:
                gauges.append(Stabilizer(gauge.only_ancilla, paulis))
        return SuperStabilizer(gauges)

    def connected(self, window: Window) -> bool:
        """Method that checks if the punctured superstabilizers are connected,
        i.e. if they share data qubits. If not, the hole will automatically be
        considered open.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.

        Output arguments:
            True or False.
        """
        gauges = [gauge for s in self.superstabilizers["initial"] for gauge in s.gauges]
        # if the hole has no gauge at all, it is not `connected``
        if not gauges:
            raise InvalidHoleFoundError(
                "Strategy dropped: Found a hole bigger than the initial patch."
            )
        # check that the hole has gauges in the initial window
        # if it has both X and Z type gauges, it's automatically `connected`
        x_gauges = [gauge for gauge in gauges if gauge.type == PauliT.X]
        z_gauges = [gauge for gauge in gauges if gauge.type == PauliT.Z]
        if any(x_gauges) and any(z_gauges):
            return True
        # if there are gauges of only one type, we need to check if
        # the data qubits of those gauges fall on the boundary which
        # explains why we miss the other type, otherwise it is not
        # `connected`. This is relevant for repurposed checks for
        # ancilla defects directly on the boundary
        data_qubits = set(q for gauge in gauges for q in gauge.data_qubits)
        return all([window.on_boundary(q) for q in data_qubits])

    def commute(self) -> bool:
        """Method that checks if the punctured superstabilizers commute.

        Output arguments:
            True or False.
        """
        return all(
            [
                stabilizers_commute(s1, s2)
                for s1 in self.superstabilizers["initial"]
                for s2 in self.superstabilizers["initial"]
            ]
        )

    def set_contours(
        self,
        data_defects: set[Pos],
        ancilla_defects: set[Pos],
        ancilla_zombies: set[Pos],
    ) -> None:
        """Method that returns the contours in the hole (everything inside the
        contours is defective/inactive.) The contours are defined only with the
        data qubits.

        Input arguments:
            `data_defects`: A set of all the data defects in the patch.
            `ancilla_defects`: A set of all the ancilla defects in the patch.
            `ancilla_zombies`: A set of all the ancilla zombies in the patch.
        """
        defects = data_defects | ancilla_defects | ancilla_zombies
        self.cycles, self.qubits_on_contours = Cycle.make_cycles(self.superstabilizers, defects)

    def is_inside(self, qubit: Pos) -> bool:
        """Method that checks if a qubit is inside the hole.

        Input arguments:
            `qubit`: Pos to inspect.

        Output arguments:
            True or False
        """
        for cycle in self.cycles:
            if cycle.polygon.contains(Point(qubit.x, qubit.y)):
                return True
        return False

    def filter_qubits_inside(
        self, window: Window, data_defects: set[Pos], ancilla_defects: set[Pos]
    ) -> None:
        """Method that filters all the qubits inside the hole per type.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
            `data_defects`: A set of all the data defects in the patch.
            `ancilla_defects`: A set of all the ancilla defects in the patch.
        """
        # all data qubits in the hole
        self.data_qubits = set(qubit for qubit in window.data_qubits if self.is_inside(qubit))
        # all ancilla qubits in the hole
        self.ancilla_qubits = set(qubit for qubit in window.ancillas if self.is_inside(qubit))
        # all data defects in the hole
        self.data_qubits_defective = self.data_qubits & data_defects
        # all ancilla defects in the hole
        self.ancilla_qubits_defective = self.ancilla_qubits & ancilla_defects
        # all zombie qubits in the hole
        self.data_qubits_zombie = self.data_qubits - self.data_qubits_defective
        # all zombie ancilla qubits
        self.ancilla_qubits_zombie = self.ancilla_qubits - self.ancilla_qubits_defective

    def clean_gauges_if_open(self, window: Window) -> None:
        """Method that returns the gauges in an open hole to be kept along the deformed
        boundary that goes around that hole.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
        """
        self.boundary_deformations = clean_gauges_if_open(
            self.cycles, window, self.superstabilizers["initial"]
        )
        assert len(self.boundary_deformations) != 0

    @classmethod
    def collect_data_to_make_holes(
        cls,
        cluster_combination: Sequence[DefectsCluster],
    ) -> tuple[list[TriangularCheck], list[SuperStabilizer], set[Pos], set[Pos], set[Pos]]:
        """Class method that collects data from the defect clusters to make holes.

        We collect the superstabilizers, defects and zombies from each cluster.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
            `cluster_combination`: A Sequence of `DefectsCluster` objects, each corresponding
            to a defects cluster in the *extended window* for which we found all possible
            strategies (i.e. Auger or efficient) corresponding to a set of superstabilizers.
            The `Hole` objects in the patch are defined by these superstabilizers.

        Output arguments:
            A tuple containing
            1) A list of `TriangularCheck` objects (efficient checks) in the patch.
            2) A list of superstabilizers.
            3) A set of data defects.
            4) A set of ancilla defects.
            5) A set of ancilla zombies.
        """
        # We try each combination (i.e. 'cluster_combination') of efficient strategies
        # for each DefectsCluster.
        # List of efficient checks
        efficient_strategy: list[TriangularCheck] = []
        # List of superstabilizers
        superstabilizers: list[SuperStabilizer] = []
        # Combine the repurposed checks and superstabilizers of all DefectClusters
        # together.
        data_defects: set[Pos] = set()
        ancilla_defects: set[Pos] = set()
        ancilla_zombies: set[Pos] = set()
        for cluster in cluster_combination:
            efficient_strategy.extend(cluster.efficient_strategy)
            superstabilizers.extend(cluster.x_superstabilizers)
            superstabilizers.extend(cluster.z_superstabilizers)
            # NB: As part of clustering, we have mapped certain link defects
            # to effective ancilla defects, and we need to collect all ancilla
            # defects including these.
            data_defects |= cluster.data_qubits_defective
            ancilla_defects |= cluster.ancillas_defective
            ancilla_zombies |= cluster.ancillas_zombie
        return efficient_strategy, superstabilizers, data_defects, ancilla_defects, ancilla_zombies

    @classmethod
    def make_holes(
        cls,
        window: Window,
        superstabilizers: list[SuperStabilizer],
        data_defects: set[Pos],
        ancilla_defects: set[Pos],
        ancilla_zombies: set[Pos],
    ) -> list[Hole]:
        """Class method that generates all holes in the patch, both closed and open.

        We use `pair_superstabilizers` to group the superstabilizers together to
        form the holes.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
            `cluster_combination`: A Sequence of `DefectsCluster` objects, each corresponding
            to a defects cluster in the *extended window* for which we found all possible
            strategies (i.e. Auger or efficient) corresponding to a set of superstabilizers.
            The `Hole` objects in the patch are defined by these superstabilizers.

        Output arguments:
            A list of `Hole` objects.
        """
        # pair superstabilizers together so we can define the holes
        x_superstabilizers = [s for s in superstabilizers if s.type == PauliT.X]
        z_superstabilizers = [s for s in superstabilizers if s.type == PauliT.Z]
        superstabilizer_pairs = pair_superstabilizers(x_superstabilizers, z_superstabilizers)
        # define the holes
        holes = [
            cls(
                window,
                x_superstabilizers,
                z_superstabilizers,
                data_defects,
                ancilla_defects,
                ancilla_zombies,
            )
            for x_superstabilizers, z_superstabilizers in superstabilizer_pairs
        ]

        return holes