"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import networkx
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from defects_module.base import (
    Pos,
    Stabilizer,
    SuperStabilizer,
)


@dataclass
class Cycle:
    """Dataclass that stores information about a cycle in a hole.
    A `Cycle` is defined by a cycle in a graph whose vertices and
    edges are determined by the superstabilizers in the hole in the
    *extended window*. Every edge correspond to a gate in one of the
    superstabilizers. The cycle must contain a defective qubit (data
    or ancilla defect) or a defective link (ancilla zombie). The cycle
    is defined by a perimeter which is a list of data qubits. The
    class also contains the gauges along the cycle in the *initial
    window*. It also stores whether or not the cycle is closed in the
    *initial window*. It also contains a `Polygon` object whose contour
    is the cycle. This polygon is used to determine the qubits that
    are contained or not inside the cycle.
    """

    qubits: list[Pos]
    """Data qubits along the cycle.

    They would ultimately be used to draw new boundaries around the
    hole if it is open.
    """
    gauges: list[Stabilizer]
    """Gauge checks along the cycle that are part of the *initial window*."""
    gauge_graph: networkx.Graph
    """Graph whose connectivity is determined by the gauge checks of
    the superstabilizers in the hole. Vertices are data qubits which
    are connected by an edge if they are part of the same gauge check.
    We only keep data qubits that are along the cycle and inside the
    *initial window*
    """
    is_closed: bool
    """Whether the cycle is closed or not in the *initial window.*"""
    polygon: Polygon
    """Polygon associated with the cycle.

    This is used to determine if a qubit is contained inside the cycle.
    """

    @staticmethod
    def _get_gauge_edges(gauge: Stabilizer) -> list[tuple[Pos, Pos]]:
        """Static method that sorts the data qubits in the gauge check
        such that they define the contour of a polygon. This is to make
        sure that we get a single cycle per gauge check, otherwise it can
        be difficult to extract the contour of a hole with a graph built
        from all the gates in gauge checks (i.e. they are many cycles).

        Input arguments:
            `gauge`: Gauge check along the cycle for which we need to
            define a contour or perimeter.

        Output arguments:
            A list of tuple of Pos. Each tuple is an edge in the cycle.
        """
        # need to have at least two points
        if len(gauge.data_qubits) < 2:
            return []
        # sort them along a circle with the ancilla at the origin
        qubits = sorted(
            gauge.data_qubits,
            key=lambda qubit: math.atan2(
                qubit.y - gauge.only_ancilla.y, qubit.x - gauge.only_ancilla.x
            ),
        )
        # return the edges along the sorted contour
        return [edge for edge in zip(qubits, qubits[1:] + [qubits[0]], strict=True)]

    @staticmethod
    def _make_graph_from_gauges(gauges: list[Stabilizer]) -> networkx.Graph:
        """Static method that returns the graph used to find the cycles
        associated with a single hole. The graph is defined from the list
        of gauge checks in the hole. For each gauge check, we define a contour
        or perimeter around the gauge check to extract a list of edges. These
        are the edges of the graph. This graph is defined in the *extended window*.

        Input arguments:
            `gauges`: A list of gauge checks.

        Output arguments:
            A networkx graph.
        """
        # define the connectivity graph of the hole used to find the contour
        # i.e. a graph composed of the data qubits found in two gauge checks
        # as nodes and the contour is defined as the largest cycle
        # create a connectivity graph for the hole based on the gauge checks connectivity
        contour_graph = networkx.from_edgelist(
            [(q1, q2) for gauge in gauges for q1, q2 in Cycle._get_gauge_edges(gauge)]
        )
        # instantiate the set of data qubits defining the contour
        contour_qubits: set[Pos] = set()
        # here we keep only the set of data qubits that are found in two
        # gauges (they are automatically on the contour of the hole)
        for gauge1, gauge2 in itertools.combinations(gauges, r=2):
            contour_qubits |= set(gauge1.data_qubits) & set(gauge2.data_qubits)
        contour_graph = contour_graph.subgraph(contour_qubits).copy()
        # return the graph
        return contour_graph

    @classmethod
    def _find_cycles(
        cls, superstabilizers: dict[str, set[SuperStabilizer]], defects: set[Pos]
    ) -> list[Cycle]:
        """Class method that finds the minimum cycle basis in a graph defined
        by the superstabilizers in the extended window and make a `Cycle` object
        for each cycle found in the hole.

        The function works as follows:
        1) For each pair of superstabilizers of opposite Pauli types in the
        extended window we create a graph with their gauges: edges are the
        gates in the gauges.
        2) For each graph we find the minimum cycle basis using networkx.
        3) For each unique cycle, we verify that the cycle encloses defects,
        by defining a Polygon object.
        4) For each gauge check, which exists in the *initial window*, we keep
        the two data qubits along the cycle.
        5) We define a graph `gauge_graph` which is defined by the gauges in
        the initial window along the cycle. We only keep the gates along the
        cycle (i.e. the graph is a line graph which can be cyclic or open).
        6) We check if `gauge_graph` is cyclic or not: if not then the cycle
        is open, otherwise it is closed in the initial window.

        Input arguments:
            `superstabilizers`: A dictionary with the sets of superstabilizers
            in both the initial and extended windows.

        Output arguments:
            A list of Cycle objects.
        """
        all_cycles: list[Cycle] = []
        for s1, s2 in itertools.combinations(superstabilizers["extended"], r=2):
            # we need two superstabilizers of opposite Pauli types
            if s2.type == s1.type:
                continue
            # create a graph from the gauges of the two superstabilizers
            subgraph = Cycle._make_graph_from_gauges(s1.gauges + s2.gauges)
            # find the cycles in the graph in the extended window
            new_cycles = networkx.minimum_cycle_basis(subgraph)
            for cycle in new_cycles:
                # we define the polygon object from the contour such as to filter
                # which qubit is contained inside the hole later
                polygon = Polygon([(pt.x, pt.y) for pt in cycle])
                # we make sure there is at least one defect inside the cycle otherwise
                # it shouldn't be around a hole and therefore we abort
                if not any([polygon.contains(Point(qubit.x, qubit.y)) for qubit in defects]):
                    continue
                # we define a list of gauges around the cycle that are in the initial window
                # and a graph whose edges are obtained from these gauges along the cycle
                # i.e. the graph is a line graph (cyclic or not)
                gauges: list[Stabilizer] = []
                gauge_graph: networkx.graph = networkx.Graph()
                # we iterate through all the gauges in the extended window
                for gauge in s1.gauges + s2.gauges:
                    # make sure the gauge is along the cycle, i.e. that it has an
                    # edge around the cycle meaning exactly two data qubits
                    # within the cycle
                    if len([q for q in gauge.data_qubits if q in cycle]) != 2:
                        continue
                    # we only keep the punctured gauges (i.e. gauges inside the initial window)
                    if gauge not in [g for s in superstabilizers["initial"] for g in s.gauges]:
                        continue
                    gauges.append(gauge)
                    gauge_graph.add_edge(*[q for q in gauge.data_qubits if q in cycle])
                # now we need to check if cycle is closed:
                # first we get the edges from the cycle itself in the *extended* window
                cycle_edges = networkx.subgraph(subgraph, nbunch=list(cycle)).edges
                # then we check that each edge in the cycle is part of the gauge checks
                # in the initial window
                is_closed = all([gauge_graph.has_edge(*edge) for edge in cycle_edges])
                # we build a Cycle and append it to the list
                all_cycles.append(Cycle(list(cycle), gauges, gauge_graph, is_closed, polygon))
        return all_cycles

    @classmethod
    def make_cycles(
        cls,
        superstabilizers: dict[str, set[SuperStabilizer]],
        defects: set[Pos],
    ) -> tuple[list[Cycle], set[Pos]]:
        """Class method that returns all the cycles in the hole
            and all the data qubits on the cycles.

        Input arguments:
            `superstabilizers`: The dictionary of superstabilizers in
            both the *initial* and *extended* windows.
            `defects`: A set of all the defective qubits. They include
             ancilla zombies that account for the link defects.

        Output arguments:
            A tuple containing
            1) A list of `Cycle` objects.
            2) A set of data qubits which live on all the contours or cycles
            in the hole.
        """
        cycles = cls._find_cycles(superstabilizers, defects)
        qubits_on_contours = {q for cycle in cycles for q in cycle.qubits}

        return cycles, qubits_on_contours
