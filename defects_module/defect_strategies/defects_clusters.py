"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx

from defects_module.base import (
    PauliT,
    Pos,
    SuperStabilizer,
)

from .make_windows import Window


@dataclass
class DefectsCluster:
    """DefectsCluster is a class that contains all the necessary information for the
    Snakes & Ladders algorithm to run on a region that includes a given cluster of
    defective qubits.

    The `graph` reflects how physical qubits around the cluster are connected. It has data
    and ancilla qubits as nodes, and its edges represent 2-qubit entangling gates between
    data and ancilla qubits.
    The graph also contains information such as the status of the qubit:
    (1) alive: the qubit is active in the code, used for both data and ancilla qubits
    (2) dead: the qubit is inactive in the code, used for both data and ancilla qubits
    (3) zombie: the qubit is disabled but could be reactivated depending on the status of its
    neighbors, only used for data qubits.
    """

    graph: networkx.Graph
    """Graph of the cluster"""
    vertical_logical: PauliT
    """Vertical logical"""
    unavailable_padding_ancillas: set[Pos]
    """Set of padding ancillas that unavailable for use."""
    x_ancillas: set[Pos] = field(default_factory=set)
    """Set of all X-type ancillas in the cluster"""
    x_ancillas_defective: set[Pos] = field(default_factory=set)
    """Set of all defective X-type ancillas in the cluster
    A defective qubit is one that is not working.
    """
    x_ancillas_active: set[Pos] = field(default_factory=set)
    """Set of all non-defective X-type ancillas in the cluster
    """
    z_ancillas: set[Pos] = field(default_factory=set)
    """Set of all Z-type ancillas in the cluster"""
    z_ancillas_defective: set[Pos] = field(default_factory=set)
    """Set of all defective Z-type ancillas in the cluster"""
    z_ancillas_active: set[Pos] = field(default_factory=set)
    """Set of all non-defective Z-type ancillas in the cluster"""
    ancillas: set[Pos] = field(default_factory=set)
    """Set of all ancillas in the cluster"""
    ancillas_defective: set[Pos] = field(default_factory=set)
    """Set of all defective ancillas in the cluster"""
    ancillas_active: set[Pos] = field(default_factory=set)
    """Set of all non-defective ancillas in the cluster"""
    data_qubits: set[Pos] = field(default_factory=set)
    """Set of all data qubits in the cluster"""
    data_qubits_defective: set[Pos] = field(default_factory=set)
    """Set of all defective data qubits in the cluster"""
    data_qubits_active: set[Pos] = field(default_factory=set)
    """Set of all non-defective data qubits in the cluster"""
    data_qubits_zombie: set[Pos] = field(default_factory=set)
    """Set of all zombie data qubits in the cluster
    A zombie qubit is not defective but currently disabled
    """
    ancillas_zombie: set[Pos] = field(default_factory=set)
    """Set of all zombie ancillas in the cluster
    A zombie ancilla is one that might be disabled because
    of a link defect.
    """
    padding_ancillas: set[Pos] = field(default_factory=set)
    """Set of padding ancillas"""
    padding_ancillas_active: set[Pos] = field(default_factory=set)
    """Set of non-defective padding ancillas"""
    padding_ancillas_defective: set[Pos] = field(default_factory=set)
    """Set of defective padding ancillas"""
    defects: set[Pos] = field(default_factory=set)
    """Set of all defective qubits in the cluster"""
    links_defective: set[tuple[Pos, Pos]] = field(default_factory=set)
    """Set of all defective links"""
    efficient_strategy: list = field(default_factory=list)
    """List of all efficient checks in the cluster"""
    x_superstabilizers: list[SuperStabilizer] = field(default_factory=list)
    """List of all X-type superstabilizers in the cluster"""
    z_superstabilizers: list[SuperStabilizer] = field(default_factory=list)
    """List of all Z-type superstabilizers in the cluster"""

    def __post_init__(self) -> None:
        """Update the fields based on the graph."""
        qubits = self.graph.nodes(data=True)

        for qubit, attrs in qubits:
            match attrs["type"]:
                case "x_ancilla":
                    self.x_ancillas.add(qubit)
                case "z_ancilla":
                    self.z_ancillas.add(qubit)
                case "data":
                    self.data_qubits.add(qubit)
                case "padding":
                    self.padding_ancillas.add(qubit)
                case _:
                    raise ValueError("Unexpected qubit status for %s" % qubit)

            match attrs["type"], attrs["status"]:
                case "x_ancilla", "dead":
                    self.x_ancillas_defective.add(qubit)
                case "z_ancilla", "dead":
                    self.z_ancillas_defective.add(qubit)
                case "x_ancilla", "alive":
                    self.x_ancillas_active.add(qubit)
                case "z_ancilla", "alive":
                    self.z_ancillas_active.add(qubit)
                case "data", "dead":
                    self.data_qubits_defective.add(qubit)
                case "data", "alive":
                    self.data_qubits_active.add(qubit)
                case "padding", "dead":
                    self.padding_ancillas_defective.add(qubit)
                case "padding", "alive":
                    self.padding_ancillas_active.add(qubit)
                case _:
                    raise ValueError("Invalid qubit attributes for %s" % qubit)

        edges = self.graph.edges(data=True)
        for node1, node2, attrs in edges:
            match attrs["defective"]:
                case False:
                    pass
                case True:
                    self.links_defective.add((node1, node2))
                case _:
                    raise ValueError("Invalid edge attributes for (%s, %s)" % (node1, node2))
        self.ancillas = self.x_ancillas | self.z_ancillas
        self.ancillas_defective = self.x_ancillas_defective | self.z_ancillas_defective
        self.ancillas_active = self.x_ancillas_active | self.z_ancillas_active
        self.defects = (
            self.ancillas_defective | self.padding_ancillas_defective | self.data_qubits_defective
        )


def generate_defects_clusters(
    window: Window,
    ancilla_defects: set[Pos],
    data_defects: set[Pos],
    link_defects: set[tuple[Pos, Pos]],
    unavailable_padding_ancillas: set[Pos],
) -> list[DefectsCluster]:
    """Function that returns a list of DefectsCluster objects for each cluster
    of defects given lists of ancilla, data, link defects in the code patch.

    Here a DefectsCluster object is a subpatch which encloses all
    defective qubits that cannot be handled by separate super-stabilizers with
    the Auger strategy. This is worst case scenario that sets an upper bound
    on the size of the super-stabilizeres. With the efficient strategy, we most
    likely end up with smaller superstabilizers for the same cluster of defects.

    To build the clusters, we start with a graph that only contains the defective qubits
    (here treat the qubits on both ends of each defective link as defects). We apply the Auger
    method to the defects to identify the surrounding qubits that need to be disabled
    or used in gauge checks, and add those qubits to the graph. The edges in the graph
    represent links between qubits. After finalizing the graph, we divide it into
    connected components so that each contains information of a cluster.

    Input arguments:
        `window`: Window object containing information about the *extended* window.
        `ancilla_defects`: A set of all ancilla defects in the patch.
        `data_defects`: A set of all data defects in the patch.
        `link_defects`: A set of all link defects in the patch.
        `unavailable_padding_ancillas`: A set of all padding ancillas
        that are either defective or are part of defective links.

    Output arguments:
        A list of `DefectsCluster` objects.
    """
    # The defective data and ancilla qubits are the first to include in the graph
    data_qubits_in_graph = data_defects.copy()
    ancillas_in_graph = ancilla_defects.copy()

    # Add qubits on both ends of the defective links to the graph
    all_window_ancillas = window.ancillas | window.padding_ancillas
    for link in link_defects:
        # ignore the defective link if it is not part of the patch
        if not set(link).issubset(window.data_qubits | all_window_ancillas):
            continue
        if link[0] in window.data_qubits:
            data, ancilla = link[0], link[1]
        else:
            data, ancilla = link[1], link[0]
        # add the data qubit on one end of the broken link to the graph
        data_qubits_in_graph.add(data)
        # add the ancilla qubit on the other end of the broken link to the graph
        ancillas_in_graph.add(ancilla)

    # Grow the clusters around the defects
    holes_graph = networkx.Graph()
    # The Auger method disables all the data qubits connected to a defective ancilla,
    # including data defects and non-defective data qubits that are turned into zombies.
    # Add these qubits to the graph.
    for ancilla in ancillas_in_graph:
        for neighbor in Pos.neighbors(ancilla):
            # if neighbor in window.data_qubits:
            data_qubits_in_graph.add(neighbor)
            holes_graph.add_edge(ancilla, neighbor)

    # For each data defect, add all its neighboring ancilla qubits
    # to the graph since we need them for the gauge checks
    for data_qubit in data_qubits_in_graph:
        for neighbor in Pos.neighbors(data_qubit):
            # if neighbor in window.ancillas | window.padding_ancillas:
            ancillas_in_graph.add(neighbor)
            holes_graph.add_edge(neighbor, data_qubit)

    # Grow the defect clusters to account for weight-1 checks
    def expand_hole(
        hole: networkx.Graph, ancillas: set[Pos]
    ) -> tuple[set[tuple[Pos, Pos]], set[Pos]]:
        """Expand the hole by including frozen data qubits and their neighboring ancillas."""
        # Define a set of edges to add to the connectivity graph because of frozen qubits
        new_edges: set[tuple[Pos, Pos]] = set()
        # Define the set of weight-1 gauge check ancillas
        dead_ancillas = {anc for anc in set(hole.nodes) & ancillas if hole.degree[anc] == 3}
        # If there are any weight-1 gauge checks, we need to zombify their data qubits
        while dead_ancillas:
            # For each weight-1 gauge check
            for dead_anc in dead_ancillas:
                # We get all the neighbors of the ancilla of the check
                for dat in Pos.neighbors(dead_anc):
                    assert dat in window.data_qubits
                    # We make sure to add all its neighbors to the graph
                    new_edges.add((dead_anc, dat))
                    hole.add_edge(dead_anc, dat)
                    # For each ancilla around those zombie data qubits
                    for anc in Pos.neighbors(dat):
                        assert anc in all_window_ancillas
                        # We make sure to add all its neighbors to the graph
                        new_edges.add((anc, dat))
                        hole.add_edge(anc, dat)
                        ancillas.add(anc)
            # We determine if there is any weight-1 gauge check left: if so, we redo the above
            dead_ancillas = {anc for anc in set(hole.nodes) & ancillas if hole.degree[anc] == 3}
        return new_edges, ancillas

    def expand_holes(
        holes_graph: networkx.Graph, ancillas_in_graph: set[Pos]
    ) -> tuple[networkx.Graph, set[Pos]]:
        """Expand all the holes in the patch by zombifying frozen qubits for the
        worst case scenario of no repurposed checks.
        """
        # define a set of all new ancillas that will be added to the connectivity
        # graph because they are around frozen qubits
        new_ancillas_in_graph: set[Pos] = set()
        # define a list of edges to add to the connectivity graph because of
        # frozen qubits
        buffer_edges: list[tuple[Pos, Pos]] = []
        # we define a hole as a connected component in the connectivity graph:
        # essentially, we connect all the gauges of both X and Z types together
        # that would be part of anti-commuting superstabilizers.
        for hole in networkx.connected_components(holes_graph):
            # for each hole, we have to determine if there are frozen qubits
            # first we identify the ancillas in the hole
            ancillas_in_hole = ancillas_in_graph & set(hole)
            # if none, we keep going
            if not any(ancillas_in_hole):
                continue
            # then we get the connectivity subgraph of the hole
            subgraph = networkx.subgraph(holes_graph, nbunch=list(hole)).copy()
            # we define a list of qubits to remove from the subgraph:
            # this is for link defects only. We will remove all functional
            # ancillas and data qubits that were added as *padding* for link
            # defects. Here we mean that if we draw a line from a weight-3
            # gauge check ancilla (on the contour of the hole) to its next
            # neighbor (another ancilla qubit inside the hole, zombie or dead)
            # that the data qubit and the ancilla in the hole are BOTH functional
            # and there is NO link defect along that line. This means that those
            # qubits are only *unnecessary padding* that will be unused when
            # repurposing.
            nodes_to_remove: set[Pos] = set()
            for node in subgraph.nodes:
                # check that the node is an ancilla
                if node not in ancillas_in_graph:
                    continue
                # check that the ancilla is on the contour
                # (here weight-1 means it connects to 1 zombie,
                # and therefore has a weight-3 check)
                if subgraph.degree[node] != 1:
                    continue
                # get its neighbouring data qubit
                data_qubit = next(subgraph.neighbors(node))
                # make sure its neighbor is not defective
                if data_qubit in data_defects:
                    continue
                # that we get the second neighbors (ancillas
                # in the hole)
                ancillas = set(subgraph.neighbors(data_qubit))
                # make sure that there is no ancilla defect
                if any(ancillas & ancilla_defects):
                    continue
                # we check if there is any link defect up to
                # the second neighbors
                links = {(anc, data_qubit) for anc in ancillas}
                links |= {(data_qubit, anc) for anc in ancillas}
                if any(links & link_defects):
                    continue
                # then those are just functional padding qubits
                # that we can remove from the current graph
                nodes_to_remove |= {data_qubit, node}
            subgraph.remove_nodes_from(nodes_to_remove)
            # now we expand the hole by getting all the frozen qubits
            # for the WORST case scenario
            new_edges, ancillas_in_hole = expand_hole(subgraph, ancillas_in_hole)
            # we update the ancillas in the connectivity graph
            new_ancillas_in_graph |= ancillas_in_hole - ancillas_in_graph
            # we add the edges to the frozen qubits
            buffer_edges += list(new_edges)
        holes_graph.add_edges_from(buffer_edges)
        return holes_graph, new_ancillas_in_graph

    # This is a recursive procedure to include the frozen qubits (from weight-1 gauge checks)
    # in the connectivity graph: we define the holes (worst-case scenario superstabilizers)
    # and determine which additional data qubits need to be zombified for them to be valid.
    while True:
        # Expand the holes
        holes_graph, new_ancillas_in_graph = expand_holes(holes_graph, ancillas_in_graph)
        # Update the ancillas in the graph
        ancillas_in_graph |= new_ancillas_in_graph
        # Break if the graph was not updated.
        if not new_ancillas_in_graph:
            break
    # Get all the holes in the graph after having added the frozen qubits
    holes = list(networkx.connected_components(holes_graph))

    # Now we have the correct zombie qubits and all the ancillas
    # that could participate in gauge checks. Build the final
    # connectivity graph from the edges and mark the defective links
    def get_edges_from_ancillas(ancillas: set[Pos]) -> list[tuple[Pos, Pos]]:
        """Return edges that connect the given ancillas to data qubits, including
        data qubits in the graph but also the outer data qubits because they are
        needed for the gauge checks.
        """
        edges: list[tuple[Pos, Pos]] = []
        # find all edges of each ancilla in the graph that connects to a data qubit
        for ancilla in ancillas:
            for neighbor in Pos.neighbors(ancilla):
                if neighbor not in window.data_qubits:
                    continue
                edges.append((ancilla, neighbor))
        return edges

    connectivity_graph = networkx.Graph()
    connectivity_graph.add_edges_from(get_edges_from_ancillas(ancillas_in_graph), defective=False)
    networkx.set_edge_attributes(
        connectivity_graph, {link: {"defective": True} for link in link_defects}
    )

    def set_node_attrs(qubits: set[Pos], qubit_type: str, status: str) -> dict[Pos, dict[str, str]]:
        """Set node attributes"""
        return {q: {"type": qubit_type, "status": status} for q in qubits}

    # Set of X ancilla defects
    x_ancilla_defects = ancilla_defects & window.x_ancillas
    # Set of Z ancilla defects
    z_ancilla_defects = ancilla_defects & window.z_ancillas
    # Set of defective padding ancillas
    padding_defects = ancilla_defects & window.padding_ancillas

    # Dict of all qubit attributes in the connectivity graph
    attrs = set_node_attrs(x_ancilla_defects, "x_ancilla", "dead")
    attrs |= set_node_attrs(window.x_ancillas - x_ancilla_defects, "x_ancilla", "alive")
    attrs |= set_node_attrs(z_ancilla_defects, "z_ancilla", "dead")
    attrs |= set_node_attrs(window.z_ancillas - z_ancilla_defects, "z_ancilla", "alive")
    attrs |= set_node_attrs(data_defects, "data", "dead")
    attrs |= set_node_attrs(window.data_qubits - data_defects, "data", "alive")
    attrs |= set_node_attrs(padding_defects, "padding", "dead")
    attrs |= set_node_attrs(window.padding_ancillas - padding_defects, "padding", "alive")

    # Add attributes to all qubits based on their type and status
    networkx.set_node_attributes(connectivity_graph, attrs)
    # we define the connected components of this connectivity graph as the
    # defect clusters. Each connected component is its own defect cluster.
    cluster_graphs = []
    for hole in holes:
        nodes = networkx.from_edgelist(get_edges_from_ancillas(hole)).nodes
        cluster_graphs.append(connectivity_graph.subgraph(nodes).copy())

    # List of defect clusters, each cluster can be sent separately to the
    # Snakes & Ladders algorithm.
    clusters = [
        DefectsCluster(graph, window.vertical_logical, unavailable_padding_ancillas)
        for graph in cluster_graphs
    ]

    return clusters