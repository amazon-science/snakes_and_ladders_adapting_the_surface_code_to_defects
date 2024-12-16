"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import networkx
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from defects_module.base import (
    Pauli,
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)


def delete_qubits_from_check(check: Stabilizer, qubits_to_remove: set[Pos]) -> Stabilizer | None:
    """Recursive function that deletes qubits from a stabilizer: it becomes recursive if the
    stabilizer is actually a superstabilizers. The function returns the new stabilizer if it
    has weight >= 2, otherwise it returns None.

    Input arguments:
        `check`: Stabilizer to be cleaned.
        `qubits_to_remove`: Set of Pos corresponding to the qubits to remove from the stabilizer.

    Ouput arguments:
        A stabilizer if at least two data qubits remain in the stabilizer, otherwise None.
    """
    # first we check if the stabilizer is actually a sueprstabilizer
    if isinstance(check, SuperStabilizer):
        # if it is a superstabilizer, we will have to clean each of its
        # gauge: we recursively call `delete_qubits_from_check` for each
        new_gauges: list[Stabilizer] = []
        for gauge in check.gauges:
            new_gauge = delete_qubits_from_check(gauge, qubits_to_remove)
            # we add the modified gauge to the list if it has at least
            # two data qubits left
            if new_gauge:
                new_gauges.append(new_gauge)
        # if there are no gauges left, then we return None
        if len(new_gauges) == 0:
            return None
        # otherwise we return the modified superstabilizer
        return SuperStabilizer(new_gauges)
    # if the stabilizer is not a superstabilizer, we then get the paulis
    # from the stabilizer that do not involve qubits to remove
    assert isinstance(check.type, PauliT)
    paulis = [
        Pauli(qubit, check.type) for qubit in check.data_qubits if qubit not in qubits_to_remove
    ]
    # if the check has less than two data qubits left, we return None
    if len(paulis) < 2:
        return None
    # otherwise we return the modified Stabilizer
    return Stabilizer(check.only_ancilla, pauli=paulis)


def remove_extrusions(
    boundaries: dict[str, list[Pos]],
    undamaged_stabilizers: set[Stabilizer],
    superstabilizers: set[SuperStabilizer],
) -> tuple[set[Stabilizer], set[SuperStabilizer]]:
    """Function that removes extrusions from the patch, i.e. all stabilizers that
    fall outside of a polygon formed by the patch boundaries.

    Input arguments:
        `boundaries`: A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as lists of Pos objects, containg information about the boundary along
        each of the four sides of the patch. The list of Pos is the list of data qubits used
        to define the boundaries.
        `undamaged_stabilizers`: A set of all the undamaged stabilizers in the patch,
        defined in the initial window.
        `superstabilizers`: A set of all the superstabilizers in the patch, also define
        in the initial window.

    Output arguments:
        A tuple containing:
        1) The set of updated undamaged stabilizers.
        2) The set of updated superstabilizers.
    """
    # get the edges along the boundaries
    edges = [
        edge
        for boundary in boundaries.values()
        for edge in zip(list(boundary), list(boundary)[1:], strict=False)
    ]
    # define a contour from the boundaries
    contour: list[Pos] = next(iter(networkx.cycle_basis(networkx.from_edgelist(edges))), [])
    # define a polygon from the contour
    polygon = Polygon([(qubit.x, qubit.y) for qubit in contour])
    # get all the data qubits in the patch
    data_qubits = set()
    for stabilizer in undamaged_stabilizers | superstabilizers:
        data_qubits |= set(stabilizer.data_qubits)
    # find all the data qubits that need to be removed
    discarded_qubits = set(
        qubit for qubit in data_qubits if not polygon.contains(Point(qubit.x, qubit.y))
    ) - set(contour)
    # find all the stabilizers that need to be modified
    modified_stabilizers = {
        s
        for s in superstabilizers | undamaged_stabilizers
        if any(set(s.data_qubits) & discarded_qubits)
    }
    # update the stabilizers after deleting the discarded data qubits
    superstabilizers -= modified_stabilizers
    undamaged_stabilizers -= modified_stabilizers
    for check in modified_stabilizers:
        # we update the stabilizer: we remove the discarded qubits
        new_check = delete_qubits_from_check(check, discarded_qubits)
        # if the new check is None, we continue
        if not new_check:
            continue
        # otherwise we check if the new check is a superstabilizer or
        # not and add the new check to the proper set
        if isinstance(new_check, SuperStabilizer):
            superstabilizers.add(new_check)
        else:
            undamaged_stabilizers.add(new_check)
    # we return the updated sets of stabilizers
    return undamaged_stabilizers, superstabilizers
