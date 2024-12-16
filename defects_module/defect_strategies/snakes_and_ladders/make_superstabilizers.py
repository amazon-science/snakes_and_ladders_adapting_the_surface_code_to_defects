"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from typing import Sequence

import networkx

from defects_module.base import (
    Pauli,
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)

from .get_efficient_strategies import (
    EfficientDefectsCluster,
)


def define_checks(cluster: EfficientDefectsCluster) -> tuple[list[Stabilizer], list[Stabilizer]]:
    """Function that list all X-type and Z-type gauge checks in the defect cluster.

    The function works as follows:
        1) For every efficient check (stored as a TriangularCheck object) we define
        a stabilizer of whose Pauli type correspond to the Pauli type of the damaged
        stabilizer that it corresponds to.
        2) For every functioning ancilla (that is not however not a zombie ancilla)
        we define a gauge check using all its active neighboring data qubits. The
        Pauli type is that of the corresponding damaged stabilizer.

    Input arguments:
        `cluster`: An EfficientDefectsCluster object which stores information about
        a strategy for a given defects cluster.

    Output arguments:
        A tuple containing
        1) A list of X-type gauge checks.
        2) A list of Z-type gauge checks.
    """
    # Initialize list of X-type checks
    x_checks: list[Stabilizer] = []
    # Initialize list of Z-type checks
    z_checks: list[Stabilizer] = []

    # Add all efficient checks, and superchecks involving padding
    for check in cluster.efficient_strategy:

        # Append efficient checks
        if check.ancilla_replaced in cluster.x_ancillas_defective | (
            cluster.x_ancillas & cluster.ancillas_zombie
        ):
            x_checks.append(check.make_stabilizer(PauliT.X))
        else:
            z_checks.append(check.make_stabilizer(PauliT.Z))

    # Add all normal and super X-type checks. The undamaged checks will be filtered out
    # when we combine the checks to form superstabilizers.
    x_checks += [
        Stabilizer(ancilla, [Pauli(qubit, PauliT.X) for qubit in cluster.graph.neighbors(ancilla)])
        for ancilla in cluster.x_ancillas_active - cluster.ancillas_zombie
    ]
    # Add all normal and gauge Z-type checks
    z_checks += [
        Stabilizer(ancilla, [Pauli(qubit, PauliT.Z) for qubit in cluster.graph.neighbors(ancilla)])
        for ancilla in cluster.z_ancillas_active - cluster.ancillas_zombie
    ]

    return x_checks, z_checks


def combine_checks(
    cluster: EfficientDefectsCluster, x_checks: list[Stabilizer], z_checks: list[Stabilizer]
) -> tuple[list[SuperStabilizer], list[SuperStabilizer]]:
    """Function that combines the X-type and Z-type gauge checks into superstabilizers based
    on their proximity to the defects.

    The function works as follows:
        1) We define two networkx graphs: one for the X-type superstabilizers and
        one for the Z-type superstabilizers.
        2) For each data qubit that is either defective or a zombie data qubit we
        extract all the neighbouring ancilla qubits (independently of their status)
        and add an edge between the data qubit and each of these ancilla in the
        graph of the Pauli type that corresponds to that of the stabilizer that
        the ancilla is associated with in the undamaged patch.
        3) For each efficient check (TriangularCheck), we add an edge between the
        replaced ancilla and replaced ancilla in the graph of the Pauli type that
        corresponds to that of the replaced ancilla in the undamaged patch. We
        also add an edge perpendicular to the efficient checks, connecting the
        weight-4 stabilizers in weight-8 supercheck, in the graph of the opposite
        Pauli type.
        4) We build a superstabilizer from each connected component in each graph.
        The Pauli type of the superstabilizer is that of the graph.

    Input arguments:
        `cluster`: An EfficientDefectsCluster object which stores information about
        a strategy for a given defects cluster.
        `x_checks`: A list of X-type gauge checks in the defects cluster.
        `z_checks`: A list of Z-type gauge checks in the defects cluster.

    Output arguments:
        A tuple containing
        1) A list of X-type superstabilizers.
        2) A list of Z-type superstabilizers.
    """
    # X check graph
    x_stabilizer_graph = networkx.Graph()
    # Z check graph
    z_stabilizer_graph = networkx.Graph()

    # Set of all disabled data qubits: data defects and zombie qubits
    data_disabled = cluster.data_qubits_defective | cluster.data_qubits_zombie

    # Add edges between gauge checks due to data defects (same as with Auger)
    for defect in data_disabled:
        # Find all direct neighbors of disabled data qubit
        anc_neighbors = set(Pos.neighbors(defect))
        for neighbor in anc_neighbors:
            if neighbor in cluster.x_ancillas:
                # Connect the X ancilla to the data defect in the X check graph
                x_stabilizer_graph.add_edge(defect, neighbor)
            elif neighbor in cluster.z_ancillas:
                # Connect the Z ancilla to the data defect in the Z check graph
                z_stabilizer_graph.add_edge(defect, neighbor)

    # Add edges due to the efficient strategy
    for check in cluster.efficient_strategy:

        # Add an edge between the defect and the repurposed ancilla
        efficient_edge = (check.ancilla, check.ancilla_replaced)

        # Find the other two ancillas in the orientation
        # perpendicular to the efficient checks which form the superchecks
        if check.data1.x == check.data2.x:
            anc1 = Pos(check.ancilla_replaced.x, check.ancilla_replaced.y - 2)
            anc2 = Pos(check.ancilla_replaced.x, check.ancilla_replaced.y + 2)
        else:
            anc1 = Pos(check.ancilla_replaced.x - 2, check.ancilla_replaced.y)
            anc2 = Pos(check.ancilla_replaced.x + 2, check.ancilla_replaced.y)
        supercheck_edge = (anc1, anc2)

        # Add edges to the detector graphs corresponding to the check type
        if check.ancilla_replaced in cluster.x_ancillas_defective | (
            cluster.x_ancillas & cluster.ancillas_zombie
        ):
            # Add an edge between the defect and the repurposed ancilla
            x_stabilizer_graph.add_edge(*efficient_edge)
            z_stabilizer_graph.add_edge(*supercheck_edge)
        else:
            # Add an edge between the defect and the repurposed ancilla
            z_stabilizer_graph.add_edge(*efficient_edge)
            x_stabilizer_graph.add_edge(*supercheck_edge)

    # Dict for gauge checks where each ancilla is assigned a single check
    checks_dict = {
        PauliT.X: dict([(check.only_ancilla, check) for check in x_checks]),
        PauliT.Z: dict([(check.only_ancilla, check) for check in z_checks]),
    }

    # Find all superstabilizers in each detector graph by taking connected components and
    # using the dicts to convert the ancillas to checks.
    def make_superstabilizer(component: Sequence[Pos], type: PauliT) -> SuperStabilizer:
        """Make a superstabilizer from a connected component in a detector graph."""
        # we get the gauges of the superstabilizer
        # NB: some zombie ancillas might be inactive and not be in checks_dict
        # Need to filter None to prevent assertion errors in `gauges`
        active_ancillas = checks_dict[type].keys()
        gauges = filter(
            None,
            [checks_dict[type].get(node, None) for node in component if node in active_ancillas],
        )
        # and build the superstabilizer
        return SuperStabilizer(list(gauges))

    # List of all X type superstabilizers
    x_superstabilizers = list(
        map(
            lambda x: make_superstabilizer(x, PauliT.X),
            networkx.connected_components(x_stabilizer_graph),
        )
    )
    # List of all Z type superstabilizers
    z_superstabilizers = list(
        map(
            lambda x: make_superstabilizer(x, PauliT.Z),
            networkx.connected_components(z_stabilizer_graph),
        )
    )

    return x_superstabilizers, z_superstabilizers
