"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass

from defects_module.base import (
    Pauli,
    PauliT,
    Pos,
    Stabilizer,
)


@dataclass
class TriangularCheck:
    """Class that stores information about an efficient check. The
    simplest case would have two opposite triangular checks around
    a single defective ancilla. More complex scenarios mix these
    triangular checks with gauge checks found in the Auger method.

    Input arguments:
        `ancilla`: The repurposed ancilla used to do the weight-2
        triangular check.
        `ancilla_replaced`: The ancilla that is being replaced, i.e.
        the ancilla associated with the damaged stabilizer that we
        are trying to reconstruct with the triangular checks.
        `data1`: The first data qubit in the triangular check.
        `data2`: The second data qubit in the triangular check.
    """

    ancilla: Pos
    """Ancilla used to measure the efficient check"""
    ancilla_replaced: Pos
    """Defective ancilla that is being replaced"""
    data1: Pos
    """First data qubit in the check"""
    data2: Pos
    """Second data qubit in the check"""

    @property
    def nodes(self) -> list[Pos]:
        """Return nodes of the check"""
        return [self.ancilla, self.data1, self.data2]

    @property
    def edges(self) -> list[tuple[Pos, Pos]]:
        """Return edges of the check"""
        return [(self.data1, self.ancilla), (self.ancilla, self.data2)]

    def make_stabilizer(self, pauli_type: PauliT) -> Stabilizer:
        """Return the gauge stabilizer for the effective check"""
        return Stabilizer(
            self.ancilla, [Pauli(self.data1, pauli_type), Pauli(self.data2, pauli_type)]
        )


@dataclass
class EfficientRepurposing:
    """Class that stores the information about all the triangular
    checks used to reconstruct a single damaged stabilizer.

    Input arguments:
        1) `replaced_ancilla`: The ancilla that is being replaced, i.e.
        the ancilla associated with the damaged stabilizer that we
        are trying to reconstruct with the triangular checks.
        2) `checks`: A list of triangular checks (efficient weight-2 checks)
        used to reconstruct the stabilizer. This list does not include
        the gauge checks found in the Auger method, only the efficient checks.
    """

    replaced_ancilla: Pos
    checks: list[TriangularCheck]
