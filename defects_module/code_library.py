"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools

from defects_module.base import (
    DistanceT,
    Logical,
    NoiseModel,
    PauliT,
    Pos,
    Stabilizer,
    X,
    Z,
)
from defects_module.translate import PhysicalCircuit

A_SCHEDULE = [Pos(1, -1), Pos(-1, -1), Pos(1, 1), Pos(-1, 1)]
B_SCHEDULE = [Pos(1, -1), Pos(1, 1), Pos(-1, -1), Pos(-1, 1)]


class RotatedSurfaceCode:
    """Logical surface code patch."""

    logicals: dict[PauliT, Logical]
    """Logicals of the code"""
    vertical_logical: PauliT
    """The Pauli type of the weight-2 checks on the top and bottom of the code.
    Also determined which type of stabilizers follow an A-type schedule."""
    stabilizer_map: dict[Pos, Stabilizer]
    """Mapping from ancilla qubits to code stabilizers."""
    noise_model: NoiseModel
    """Map from pairs (gate, qubit) to a list of noise additions."""

    def __init__(
        self,
        # Horizontal and vertical dimensions of the code.
        distance: DistanceT,
        # Grid position of the southwest corner of the code.
        sw_offset: Pos = Pos(1, 1),
        # Pauli type of the bulk stabilizer in the southwest corner of the code.
        sw_type: PauliT = PauliT.Z,
        # Pauli type of the weight-2 stabilizers along the top and bottom of the code.
        vertical_logical: PauliT = PauliT.X,
    ):
        self._init_attributes(distance, sw_offset, sw_type, vertical_logical)
        self.stabilizer_map = self._make_stabilizers()
        self._stabilizers = set(self.stabilizer_map.values())
        self.logicals = self._make_logicals()

    def _init_attributes(
        self,
        distance: DistanceT,
        sw_offset: Pos = Pos(1, 1),
        sw_type: PauliT = PauliT.Z,
        vertical_logical: PauliT = PauliT.X,
    ):
        self.distance = (distance, distance) if isinstance(distance, int) else distance
        """Distance of the code"""

        origin = sw_offset - Pos(1, 1)
        self.extent = (
            origin.x,
            origin.x + 2 * self.distance[0],
            origin.y,
            origin.y + 2 * self.distance[1],
        )
        """Spatial extent of the code"""

        self.sw_offset = sw_offset
        """Southwest stabilizer offset"""
        self.sw_type = sw_type
        """Southwest check type"""
        self.vertical_logical = vertical_logical
        """Vertical logical"""

    def _make_stabilizers(self) -> dict[Pos, Stabilizer]:
        """Makes the stabilizer map"""
        stabilizer_map = {}
        xmin, xmax, ymin, ymax = self.extent
        # Iterate through all ancilla positions on the grid.
        for x, y in itertools.product(range(xmin, xmax + 1, 2), range(ymin, ymax + 1, 2)):
            pauli_type = self._get_stabilizer_type(x, y)

            # Skip ancillas which are unused
            if x in (xmin, xmax):
                if y in (ymin, ymax):
                    continue
                if pauli_type == self.vertical_logical:
                    continue
            if y in (ymin, ymax):
                if pauli_type != self.vertical_logical:
                    continue

            # Build the stabilizer
            anc = Pos(x, y)
            pauli_op = X if pauli_type == PauliT.X else Z
            anc_neighbors = Pos.neighbors(anc)
            data_qubits = [q for q in anc_neighbors if xmin <= q.x <= xmax and ymin <= q.y <= ymax]
            stabilizer = Stabilizer(anc, [pauli_op(q) for q in data_qubits])

            stabilizer_map[anc] = stabilizer

        return stabilizer_map

    def _make_logicals(self) -> dict[PauliT, Logical]:
        """Makes the logicals."""
        xmin, xmax, ymin, ymax = self.extent
        P = self.vertical_logical
        Q = PauliT.X if self.vertical_logical == PauliT.Z else PauliT.Z
        return {
            P: [Pos(xmin + 1, y) for y in range(ymin + 1, ymax, 2)],
            Q: [Pos(x, ymin + 1) for x in range(xmin + 1, xmax, 2)],
        }

    def _get_stabilizer_type(self, x: int, y: int) -> PauliT:
        """Get the stabilizer type at position (x, y)"""
        xmin, _, ymin, _ = self.extent
        if not (x - xmin) % 2 == (y - ymin) % 2 == 0:
            raise ValueError("Not a valid stabilizer position.")
        other_type = PauliT.X if self.sw_type == PauliT.Z else PauliT.Z
        is_sw_type = (x - xmin + y - ymin) % 4 == 0
        return self.sw_type if is_sw_type else other_type

    @property
    def stabilizers(self) -> list[Stabilizer]:
        """All stabilizers of the patch, including super base and super stabilizers.

        We return this in sorted order to ensure programs have deterministic ordering.
        """
        # XXX: Would it be better to make stabilizers, data qubits, and
        #      ancilla qubits come natively as sets, rather than lists?
        return sorted(self._stabilizers, key=lambda stab: stab.ancilla[0])

    @property
    def data_qubits(self):
        """Data qubits"""
        data_qubits = set()
        for stabilizer in self.stabilizers:
            data_qubits.update(set(stabilizer.data_qubits))

        return sorted(data_qubits)

    @property
    def ancilla_qubits(self):
        """Ancilla qubits"""
        ancilla_qubits = []
        for stabilizer in self.stabilizers:
            ancilla_qubits.extend(list(stabilizer.ancilla))
        return sorted(ancilla_qubits)

    def _check(self, circ: PhysicalCircuit, stabilizers: list[Stabilizer] | None = None):
        """Performs a check round, consisting of 8 timesteps."""
        if not stabilizers:
            stabilizers = self.stabilizers
        # XXX: Should we make class properties ancilla_qubits_X and ancilla_qubits_Z?
        ancillas_X = [anc for stab in stabilizers for anc in stab.ancilla if stab.type == PauliT.X]
        ancillas_Z = [anc for stab in stabilizers for anc in stab.ancilla if stab.type == PauliT.Z]

        # Reset
        circ.add_gate("R", ancillas_X + ancillas_Z)
        circ.add_gate("I", self.data_qubits)

        # Hadamards
        circ.add_gate("H", ancillas_X)
        circ.add_gate("I", ancillas_Z)
        circ.add_gate("I", self.data_qubits)

        # Entangling gates
        for timestep in range(4):
            idling_qubits = set(self.data_qubits) | set(ancillas_X + ancillas_Z)

            # Add a gate for each active stabilizer
            for stab in stabilizers:
                schedule = A_SCHEDULE if stab.type == self.vertical_logical else B_SCHEDULE
                # XXX: Assume here that there's only one ancilla.
                (anc,) = stab.ancilla
                data = anc + schedule[timestep]
                if data in stab.data_qubits:
                    idling_qubits -= {anc, data}
                    targets = (anc, data) if stab.type == PauliT.X else (data, anc)
                    circ.add_gate("CX", targets)

            # Add idling instructions to idling qubits
            circ.add_gate("I", list(idling_qubits))

        # Hadamards
        circ.add_gate("H", ancillas_X)
        circ.add_gate("I", ancillas_Z)
        circ.add_gate("I", self.data_qubits)

        # Measure
        circ.add_gate("M", ancillas_X + ancillas_Z)
        circ.add_gate("I", self.data_qubits)

    def initialize(self, circ: PhysicalCircuit, paulit: PauliT):
        """Initializes the logical patch in the given Pauli basis, then performs
        a single round of checks.
        """
        # Reset the qubits
        circ.add_gate("R", self.data_qubits)

        # Rotate basis
        match paulit:
            case PauliT.Z:
                pass
            case PauliT.X:
                circ.add_gate("H", self.data_qubits)
            case _:
                raise NotImplementedError

        # Do checks
        self._check(circ)

        # Add detectors at deterministic locations
        for anc, stab in self.stabilizer_map.items():
            if stab.type != paulit:
                continue
            circ.add_detector([(anc, -1)])

    def measure(self, circ: PhysicalCircuit, paulit: PauliT):
        """Measures the logical patch in the given Pauli basis."""
        # Rotate basis.
        match paulit:
            case PauliT.Z:
                pass
            case PauliT.X:
                circ.add_gate("H", self.data_qubits)
            case _:
                raise NotImplementedError

        # Measure qubits
        circ.add_gate("M", self.data_qubits)

        # Add detectors at stabilizers of that type
        for anc, stab in self.stabilizer_map.items():
            if stab.type != paulit:
                continue
            circ.add_detector([(anc, -1)] + [(q, -1) for q in stab.data_qubits])

        # Add observable.
        circ.add_observable(self.logicals[paulit])

    def check(self, circ: PhysicalCircuit):
        """Performs one check round."""
        # Do checks
        self._check(circ=circ)

        # Add detectors
        for anc in self.ancilla_qubits:
            circ.add_detector([(anc, -1), (anc, -2)])
