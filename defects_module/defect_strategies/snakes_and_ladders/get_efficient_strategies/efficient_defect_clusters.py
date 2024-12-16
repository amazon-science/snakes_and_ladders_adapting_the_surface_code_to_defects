"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Sequence

import networkx

from defects_module.base import (
    PauliT,
    Pos,
)
from defects_module.defect_strategies.defects_clusters import DefectsCluster

from .define_efficient_checks import EfficientRepurposing, TriangularCheck


@dataclass
class EfficientDefectsCluster(DefectsCluster):
    """Class that contains information about a DefectsCluster on which we applied
    efficient checks.

    Input arguments:
        `cluster`: A DefectsCluster object storing information about the defects
        and the zombies in the defect cluster.
    """

    def __init__(self, cluster: DefectsCluster) -> None:
        self.__dict__.update(cluster.__dict__)
        self.efficient_strategy: list[TriangularCheck] = []
        self.efficient_checks_graph: networkx.classes.graph.Graph | None = None
        x = [anc.x for anc in self.ancillas]
        y = [anc.y for anc in self.ancillas]
        self.xlims = (min(x) - 2, max(x) + 2)
        self.ylims = (min(y) - 2, max(y) + 2)
        self.distance_loss: tuple[int, int] | None = None

    def __equal_strategies__(self, cluster: EfficientDefectsCluster) -> bool:
        """Method that checks if two clusters have the same strategies.

        Input arguments:
            `cluster`: The other EfficientDefectsCluster to compare against.

        Output arguments:
            True if they have the same efficient strategy, else False.
        """
        checks1, checks2 = self.efficient_strategy.copy(), cluster.efficient_strategy.copy()
        # check if we have the same number of checks
        if len(checks1) != len(checks2):
            return False
        # check we have the same sets active data qubits
        if self.data_qubits_active != cluster.data_qubits_active:
            return False
        # next we check if every check in the first list is in the
        # second list (if not, then they are not equivalent)
        for check1 in checks1:
            # we define a bool flag to determine if in the list
            in_list = False
            for check2 in checks2:
                # we check if two checks are the same by
                # comparing their nodes
                if set(check1.nodes) == set(check2.nodes):
                    in_list = True
                    break
            # if not in the list then the checks are not equivalent
            if not in_list:
                return False
            # if in the list then we remove the check from the
            # second list to speed up the rest of the search
            checks2.remove(check2)
        # if we reach this point, the checks are equivalent
        return True

    def __equal_superstabilizers__(self, cluster: EfficientDefectsCluster) -> bool:
        """Method that checks if two clusters have the same superstabilizers.

        Input arguments:
            `cluster`: The other EfficientDefectsCluster to compare against.

        Output arguments:
            True if they have the same superstabilizers, else False.
        """
        s1 = self.x_superstabilizers + self.z_superstabilizers
        s2 = cluster.x_superstabilizers + cluster.z_superstabilizers
        # we check if we have the same superstabilizers
        return set(s1) == set(s2)

    @staticmethod
    def _get_sentenced_qubits_from_missing_checks(repurposing: EfficientRepurposing) -> set[Pos]:
        """Helper method that gets the data qubits that need to be disabled around
        an ancilla that is to be partially or totally replaced.

        Input arguments:
            `repurposing`: EfficientRepurposing object that contains the efficient
            checks and the replaced ancilla defect / zombie.

        Output arguments:
            A set of data qubits that need to be zombified.
        """
        anc, checks = repurposing.replaced_ancilla, repurposing.checks
        neighbors = {anc + s for s in [Pos(-1, -1), Pos(-1, 1), Pos(1, -1), Pos(1, 1)]}
        check_nodes = {q for check in checks for q in check.nodes}
        return neighbors - check_nodes

    def zombify_data_qubit_from_missing_checks(
        self,
        checks_comb: Sequence[EfficientRepurposing] | None,
    ) -> None:
        """Method that transforms data qubits into zombies. For each given combination
        of efficient repurposing, validate the efficient checks. For the efficient checks
        that don't work due to defects, remove them and set appropriate strings of qubits
        to zombies.

        Input arguments:
            `checks_comb`: A list of EfficientRepurposingObjects, which contains a list of
            all the replaced ancillas and the efficient checks around them.
        """
        # in the case of Auger, we will create a list of empty repurposings
        if checks_comb is None:
            # determine all ancillas that can be supposedly replaced
            replaced_ancillas = self.ancillas_defective | self.ancillas_zombie
            checks_comb = [EfficientRepurposing(anc, []) for anc in replaced_ancillas]

        # Define sentenced qubits as all zombie qubits that are neighbors
        # to an ancilla defect or zombie but not in the checks around it
        sentenced_qubits = {
            q
            for repurposing in checks_comb
            for q in self._get_sentenced_qubits_from_missing_checks(repurposing)
        }

        # Get all the triangular checks
        checks = [check for repurposing in checks_comb for check in repurposing.checks]

        # Create a graph for the efficient checks
        efficient_checks_graph = networkx.Graph()
        # Add the qubits and links in all the possible efficient checks to the graph
        edges = [pair for check in checks for pair in check.edges]
        efficient_checks_graph.add_edges_from(edges)

        # diagonal weight-2 checks don't count as repurposed
        ancilla_zombies = set(
            check.ancilla for check in checks if check.ancilla == check.ancilla_replaced
        )
        # List all repurposed ancillas
        repurposed_ancillas = set(check.ancilla for check in checks) - ancilla_zombies

        # Strings (connected components) in the graph that involve any defect
        def is_defective(snake: Sequence[Pos]) -> bool:
            # make sure there is no zombie or defect in the checks
            if not set(snake).isdisjoint(sentenced_qubits):
                return True
            # Confirm that all repurposed ancillas (ancilla in any check) are
            # repurposed only once by checking that the degree of each ancilla in
            # the graph is two. The strategy is only valid if that is true. HOWEVER,
            # because of how we set up the algorithm, we need to allow zombie
            # ancillas to be repurposed once or twice, since they can be still active
            # in their native check. If not valid, return an empty list for the efficient
            # strategies.
            # condition on all active ancillas that are NOT zombies: need to have degree 2
            subgraph = efficient_checks_graph.subgraph(list(snake))
            cond1 = all(
                subgraph.degree[node] == 2 for node in subgraph.nodes if node in repurposed_ancillas
            )
            # condition on the zombie ancillas: degree 2 or 3
            cond2 = all(
                2 <= subgraph.degree[node] <= 3
                for node in subgraph.nodes
                if node in ancilla_zombies
            )
            if not (cond1 and cond2):
                return True
            # check that we don't have two overlapping
            # repurposed checks
            for node in subgraph.nodes:
                if node not in ancilla_zombies:
                    continue
                if subgraph.degree[node] != 3:
                    continue
                count = 0
                for check in checks:
                    if check.ancilla == node:
                        count += 1
                if count > 2:
                    return True
            return False

        snakes = [
            snake
            for snake in networkx.connected_components(efficient_checks_graph)
            if is_defective(snake)
        ]
        # All qubits marked as sentenced because they are part of strings containing
        # defects.
        sentenced_qubits |= set(itertools.chain(*snakes))

        # List of all efficient checks that are free of sentenced qubits
        efficient_strategy = [
            check for check in checks if set(check.nodes).isdisjoint(sentenced_qubits)
        ]
        # Remove all the sentenced qubits from the graph
        efficient_checks_graph.remove_nodes_from(sentenced_qubits)

        # Store the efficient checks
        self.efficient_strategy = efficient_strategy
        # Store the efficient checks graph
        self.efficient_checks_graph = efficient_checks_graph
        # Update the zombie data qubits
        self.data_qubits_zombie = sentenced_qubits & self.data_qubits
        # Update the status of data qubit in the cluster
        self.data_qubits_active = self.data_qubits - (
            self.data_qubits_defective | self.data_qubits_zombie
        )

    def _get_frozen_qubits(self) -> set[Pos]:
        """Helper method that extracts the frozen qubits in the cluster. By
        frozen, we mean data qubits that result in weight-1 checks.

        Output arguments:
            A set of qubits that are frozen.
        """
        frozen_qubits: set[Pos] = set()
        for anc in self.ancillas_active:
            # check for weight-1 gauge checks which we remove
            if self.graph.degree[anc] == 1:
                # the data qubits involved in those are 'frozen' and
                # need to be killed
                frozen_qubits.add(next(self.graph.neighbors(anc)))
        return frozen_qubits

    def _clean_invalid_efficient_checks(self, qubits_to_remove: set[Pos]) -> set[Pos]:
        """Helper method that removes efficient checks with dead qubits.

        Input arguments:
            `qubits_to_remove`: A set of qubits that are zombified or defective.
            We use them to remove invalid efficient checks. If the efficient
            check is removed, all its data qubits are zombified.

        Output arguments:
            A set of data qubits that are zombified.
        """
        # Strings (connected components) in the graph that involve any defect
        def is_defective(snake: Sequence[Pos]) -> bool:
            # make sure there is no zombie or defect in the checks
            if not set(snake).isdisjoint(qubits_to_remove):
                return True
            return False

        snakes = [
            snake
            for snake in networkx.connected_components(self.efficient_checks_graph)
            if is_defective(snake)
        ]
        # All qubits marked as sentenced because they are part of strings containing
        # defects.
        sentenced_qubits = set(itertools.chain(*snakes))

        # List of all efficient checks that are free of sentenced qubits
        self.efficient_strategy = [
            check
            for check in self.efficient_strategy
            if set(check.nodes).isdisjoint(sentenced_qubits)
        ]

        # Return the data qubits that are zombified
        return sentenced_qubits & self.data_qubits

    def _remove_sentenced_qubits_from_graph(self, qubits_to_remove: set[Pos]) -> None:
        """Helper method that removes the sentenced qubits from the cluster graph.

        Input arguments:
            `qubits_to_remove`: A set of frozen/zombie qubits to remove from the graph.
        """
        # Remove all edges to the remaining zombie qubits. This is because they are now
        # effectively defective.
        self.data_qubits_zombie |= qubits_to_remove - self.data_qubits_defective
        # Remove all edges to the remaining zombie qubits. This is because they are now
        # effectively defective.
        edges_to_remove = [
            (neighbor, defect)
            for defect in qubits_to_remove
            for neighbor in self.graph.neighbors(defect)
        ]
        self.graph.remove_edges_from(edges_to_remove)

    def _zombify_decoupled_ancillas(self) -> None:
        """Helper method that sets the zombie ancillas of the defects cluster."""
        for anc in self.ancillas_active:
            if self.graph.degree[anc] == 0:
                self.ancillas_zombie.add(anc)

    def zombify_frozen_qubits(self) -> None:
        """Method that removes zombie data qubits and recursively remove frozen qubits."""
        sentenced_qubits = self.data_qubits_zombie | self.data_qubits_defective
        while sentenced_qubits:
            self._remove_sentenced_qubits_from_graph(sentenced_qubits)
            sentenced_qubits = self._get_frozen_qubits()
            sentenced_qubits |= self._clean_invalid_efficient_checks(sentenced_qubits)
        # Update the status of data qubits in the cluster
        self.data_qubits_active = self.data_qubits - (
            self.data_qubits_defective | self.data_qubits_zombie
        )
        self._zombify_decoupled_ancillas()

    def _compute_distance_loss_per_pauli_type(
        self, pauli_type: PauliT, boundary_nodes: tuple[str, str], lims: tuple[int, int], coord: str
    ) -> tuple[int, list[list[list[Pos]]]]:
        """Helper method that computes the distance for a given Pauli type, a set
        of virtual boundary nodes, the limits of the region of the graph and the coordinate
        (orientation) of the graph.

        The function works as follows:
            1) For each efficient check we either add an edge between the replaced ancilla
            and the repurposed ancilla if the Pauli type of the replaced ancilla is that of
            the graph, or an edge perpendicular to the efficient check between the two
            neighbouring ancillas found in the supercheck.
            2) For each data defect or data zombie we add an edge of weight 0 to every of its
            neighbouring ancilla qubits.
            3) For every pair of ancillas in the graph of the Pauli type `pauli_type` we add an
            edge between them with a weight corresponding to the Chebyshev distance between them.
            4) For every ancilla in the graph for the Pauli_type `pauli_type`, we connect it to
            the boundary nodes still using the Chebyshev distance.
            5) We use the dijkstra algorithm to compute the shortest distance between the
            boundary nodes.

        Input arguments:
            `pauli_type`: Pauli type of the logicals for whcih we want to compute the shortest
            distance.
            `boundary_nodes`: Virtual boundary nodes.
            `lims`: Limits of the cluster (i.e. size of the cluster) along the axis of the
            logicals to compute.
            `coord`: Coordinate 'x' (horizontal) or 'y' (vertical) corresponding to the axis
            of the logicals to compute.

        Output arguments:
            An integer corresponding to the distance loss.
            The set of data qubits in each gauge, in each superstabilizer to check the validity
            of the solution.
        """
        # list of data qubits in each superstabilizer: they are divided as sets, one for each
        # gauge in the superstabilizer
        data_qubits: list[list[list[Pos]]] = []
        # ancillas in the detector graph
        ancillas = self.x_ancillas if pauli_type == PauliT.X else self.z_ancillas
        # detector graph
        graph = networkx.Graph()
        # add edges due to the efficient checks
        for check in self.efficient_strategy:
            if check.ancilla_replaced in ancillas:
                edge = (check.ancilla, check.ancilla_replaced)
            else:
                if check.data1.x == check.data2.x:
                    anc1 = Pos(check.ancilla_replaced.x, check.ancilla_replaced.y - 2)
                    anc2 = Pos(check.ancilla_replaced.x, check.ancilla_replaced.y + 2)
                else:
                    anc1 = Pos(check.ancilla_replaced.x - 2, check.ancilla_replaced.y)
                    anc2 = Pos(check.ancilla_replaced.x + 2, check.ancilla_replaced.y)
                edge = (anc1, anc2)
            graph.add_edge(*edge, weight=0.0)
        # add edges due to data qubit defects
        for q in self.data_qubits_zombie | self.data_qubits_defective:
            neighbors = [
                q + shift
                for shift in [Pos(-1, -1), Pos(-1, 1), Pos(1, -1), Pos(1, 1)]
                if q + shift in ancillas
            ]
            for edge in itertools.combinations(neighbors, r=2):
                graph.add_edge(*edge, weight=0.0)

        # Let's make sure that the superstabilizers are valid:
        # make sure that each data qubit in the superstabilizer
        # is acted on by only one single gauge check in the
        # superstabilizer. This is to avoid the case where the
        # superstabilizer folds onto itself.
        for component in networkx.connected_components(graph):
            # for each connected component in the graph (here
            # containing only weight-0 edges between ancillas
            # neighboring defects), which corresponds to its
            # own superstabilizer, we get the data qubits in
            # all the gauges of the superstabilizer
            data_qubits_in_super: list[list[Pos]] = []
            for anc in set(component):
                # we use the ancillas that are not repurposed
                if anc not in ancillas:
                    continue
                # and keep all the data qubits that are their
                # neighbors and ALIVE: they are the data qubits
                # of a single gauge
                data_qubits_in_super.append(
                    [dat for dat in self.graph.neighbors(anc) if dat in self.data_qubits_active]
                )
            # Then we check that each data qubit in the
            # superstabilizer is only acted on once.
            dat = [d for g in data_qubits_in_super for d in g]
            for d in set(dat):
                # if the data qubit is acted on more than once,
                # then we return -1 for an invalid strategy.
                if dat.count(d) > 1:
                    return -1, []
            # we append the data qubits of the superstabilizer
            data_qubits.append(data_qubits_in_super)

        # add the edges between the detectors that are not part of the same holes
        for node1, node2 in itertools.combinations(graph.nodes, r=2):
            if (node1, node2) not in graph.edges and node1 in ancillas and node2 in ancillas:
                graph.add_edge(
                    node1, node2, weight=max(abs(node2.x - node1.x), abs(node2.y - node1.y)) / 2
                )
        # add edges to the boundary nodes
        for node in ancillas:
            # We connect all nodes to the boundary nodes with
            # weights corresponding to the shortest distance to them.
            graph.add_edge(boundary_nodes[0], node, weight=(getattr(node, coord) - lims[0]) / 2)
            graph.add_edge(boundary_nodes[1], node, weight=(lims[1] - getattr(node, coord)) / 2)

        distance = networkx.dijkstra_path_length(
            graph, boundary_nodes[0], boundary_nodes[1], weight="weight"
        )
        # compute the distance loss
        return int((lims[1] - lims[0]) / 2 - distance), data_qubits

    def compute_distance_loss(self) -> tuple[int, int]:
        """Method that estimates the distance loss of a cluster.

        This function calls the `_compute_distance_loss_per_pauli_type` method above
        for X and Z pauli types (i.e. for both horizontal and vertical logicals).

        Output arguments:
            A tuple of integers corresponding to the distance loss in both directions.
        """
        # Determine graphs for horizontal and vertical distances
        if self.vertical_logical == PauliT.X:
            hor_dist_loss, hor_dat = self._compute_distance_loss_per_pauli_type(
                PauliT.X, ("L", "R"), self.xlims, "x"
            )
            ver_dist_loss, ver_dat = self._compute_distance_loss_per_pauli_type(
                PauliT.Z, ("B", "T"), self.ylims, "y"
            )
        else:
            hor_dist_loss, hor_dat = self._compute_distance_loss_per_pauli_type(
                PauliT.Z, ("L", "R"), self.xlims, "x"
            )
            ver_dist_loss, ver_dat = self._compute_distance_loss_per_pauli_type(
                PauliT.X, ("B", "T"), self.ylims, "y"
            )
        # Now we validate that the superstabilizers commute: we check that each
        # gauge of superstabilizer A commutes with superstabilizer B
        for h_super in hor_dat:
            for v_super in ver_dat:
                # check the horizontal type gauges
                for h_gauge in h_super:
                    if len(set(h_gauge) & {q for v_gauge in v_super for q in v_gauge}) % 2:
                        # if it doesn't commute, the solution is invalid
                        return (-1, -1)
                # check the vertical type gauges
                for v_gauge in v_super:
                    if len(set(v_gauge) & {q for h_gauge in h_super for q in h_gauge}) % 2:
                        # if it doesn't commute, the solution is invalid
                        return (-1, -1)
        # return the distance loss in both directions
        return hor_dist_loss, ver_dist_loss
