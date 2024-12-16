"""
Provides the base classes and methods necessary to define surface code patches,
use them to write memory experiments in Stim, and then simulate + decode to
obtain logical error rates.

Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence, Union

###
### Base classes and types.
###


class Pos(int):
    """Pos class that subclasses int to manage positions in circuit writers,
    using a a pairing function to achieve a unique integer mapping.

    Warning: Negative coordinates are handled in a restrictive and special way, and their use in
    contexts where integers are expected is discouraged. If either x or y is negative, the
    coordinate is mapped to a negative integer. The class also imposes the following constraints
    and mappings:

    ‣ Both x and y must be greater than -1000.
    ‣ NW quadrant (x < 0 and y >= 0) maps to (0, -1e6].
    ‣ SE quadrant (x >= 0 and y < 0) maps to (-1e6, -2e6].
    ‣ SW quadrant (x < 0 and y < 0) maps to (-2e6, -3e6].

    Alternative strategies for negative coordinates, such as mapping all integers (positive and
    negative) to a non-negative sequence, were considered. This method involves transforming each
    number x to 2x if x is non-negative, and -2x-1 if x is negative, which achieves only 50% packing
    efficiency. The lower efficiency would require the allocation of more qubits in simulators,
    resulting in a performance degradation. Thus, this approach was not adopted.
    """

    MAX_NEG = -1000
    """ Maximum allowed negative coordinate value """
    NEG_QUAD_SHIFT = 1_000_000
    """ Shift amount for negative quadrants to prevent overlap. This number is based on the fact
    that Pos(999, 999) => 999999
    """

    def __new__(cls, x: int, y: int, z: int = 0):

        # The largest negative coordinate is constrainted
        assert (
            x > Pos.MAX_NEG and y > Pos.MAX_NEG
        ), f"{x =} and {y =} should must be greater than {Pos.MAX_NEG}"

        def elegant_pair(x, y):
            if x >= y:
                return x**2 + x + y
            return y**2 + x

        # calculate the integer value assuming |x| and |y|
        value = elegant_pair(abs(x), abs(y))

        # NW quadrant occupies (0, 1e6]
        if x < 0 and y >= 0:
            value = -value
        # SE quadrant occupies (-1e6, -2e6]
        elif y < 0 and x >= 0:
            value = -value - Pos.NEG_QUAD_SHIFT
        # SW quadrant occupies (-2e6, -3e6]
        elif x < 0 and y < 0:
            value = -value - 2 * Pos.NEG_QUAD_SHIFT

        object = int.__new__(cls, value)
        return object

    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y

    def __sub__(self, other: int | Pos) -> Pos:
        if isinstance(other, Pos):
            return Pos(self.x - other.x, self.y - other.y)
        return Pos(self.x - other, self.y)

    def __add__(self, other: int | Pos) -> Pos:
        if isinstance(other, Pos):
            return Pos(self.x + other.x, self.y + other.y)
        return Pos(self.x - other, self.y)

    @staticmethod
    def neighbors(pos: Pos) -> list[Pos]:
        return [
            Pos(pos.x - 1, pos.y + 1),
            Pos(pos.x + 1, pos.y + 1),
            Pos(pos.x - 1, pos.y - 1),
            Pos(pos.x + 1, pos.y - 1),
        ]

    def __repr__(self) -> str:
        return f"{self.x, self.y}"

    def __reduce__(self):
        """Custom pickling of int subclass to enable multiprocessing using Dask."""
        return (Pos, (self.x, self.y))

    def __setstate__(self, state):
        """Custom unpickling of int subclass to enable multiprocessing using Dask."""
        self.x, self.y = state


class PauliT(Enum):
    """A Pauli type with ± versions encoded as a 3-bit value"""

    Z = 0b010
    MinusZ = 0b110
    X = 0b001
    MinusX = 0b101
    Id = 0b000

    def __str__(self):
        return f"{self.name}"

    def anticommute(self, other_pauli: PauliT) -> bool:
        """Returns False if the pauli's commute, and True if they anticommute.
        Done by checking whether the last two bits are identical.
        """
        return ((self.value ^ other_pauli.value) % 4) != 0


@dataclass(frozen=True)
class Pauli:
    qubit: int | Pos
    type: PauliT

    def __repr__(self) -> str:
        return f"{self.type.name}{self.qubit}"


def _pauli(*args: int, type: PauliT):
    if len(args) == 1:
        return Pauli(args[0], type)
    return Pauli(Pos(args[0], args[1]), type)


def Z(*args: int) -> Pauli:
    return _pauli(*args, type=PauliT.Z)


def X(*args: int) -> Pauli:
    return _pauli(*args, type=PauliT.X)


class Stabilizer:
    """Stabilizer class."""

    def __init__(self, ancilla: Pos | Sequence[Pos], pauli: Sequence[Pauli]) -> None:
        self._ancilla = tuple(ancilla) if isinstance(ancilla, Sequence) else (ancilla,)
        self._pauli = tuple(pauli)

    @property
    def ancilla(self):
        return self._ancilla

    @property
    def only_ancilla(self) -> Pos:
        """If there is only a single ancilla position, return it. Otherwise, raise an error."""
        if len(self.ancilla) != 1:
            raise ValueError(f"Expected a single ancilla, but have {len(self.ancilla)}")
        return self.ancilla[0]

    @property
    def pauli(self):
        return self._pauli

    @property
    def pauli_str(self):
        return "".join([str(_pauli.type) if _pauli else "I" for _pauli in self.pauli])

    @property
    def data_qubits(self):
        return [pauli.qubit for pauli in self.pauli]

    @property
    def type(self) -> PauliT | None:
        """Return the Pauli type if all data qubit Paulis are of the same type, else None."""
        pauli_word = set(pauli.type for pauli in self.pauli)
        if len(pauli_word) == 1:
            return pauli_word.pop()
        return None

    @property
    def weight(self) -> int:
        """Number of non-identity Pauli terms"""
        return len(self.data_qubits)

    def __repr__(self) -> str:
        repr = ""
        for pauli in self.pauli:
            repr += str(pauli) if pauli else "I(∅)"
        return repr

    def __hash__(self) -> int:
        return hash((*self.ancilla, *sorted(self.pauli, key=lambda pauli: pauli.qubit)))

    def __eq__(self, other: object) -> bool:
        return hash(self) == hash(other)


class SuperStabilizer(Stabilizer):
    """Super stabilizer class holding a list of check references.

    Note that superstabilizers are a *generalization* of stabilizers -- every
    stabilizer is a superstabilizer with a single gauge. In building detectors
    for surface codes with defects, it will be convenient to call undamaged
    stabilizers whose ancillas are serving dual-purpose as superstabilizers.
    """

    def __init__(self, gauges: list[Stabilizer]) -> None:
        paulis = set()
        ancillas = set()
        for gauge in gauges:
            paulis |= set(gauge.pauli)
            ancillas |= set(gauge.ancilla)

        sorted_paulis = sorted(paulis, key=lambda pauli: pauli.qubit if pauli else 0)
        sorted_ancillas = sorted(ancillas)

        super().__init__(tuple(sorted_ancillas), tuple(sorted_paulis))
        self.gauges = gauges

    @staticmethod
    def decompose(stabilizers: Sequence[Stabilizer]) -> list[Stabilizer]:
        """Decompose any superstabilizers in the incoming list into gauges."""
        stabilizers_and_gauges = []
        for stabilizer in stabilizers:
            if isinstance(stabilizer, SuperStabilizer):
                stabilizers_and_gauges += stabilizer.gauges
            else:
                stabilizers_and_gauges.append(stabilizer)

        return stabilizers_and_gauges


GateT = str  # e.g. H, CZ, CX, ...
QubitT = Pos
NoiseT = tuple[str, float]  # e.g. DEPOLARIZE(0.1), ...
NoiseModel = dict[GateT, list[NoiseT]]  # TODO Make this qubit-specific.
DistanceT = Union[int, tuple[int, int]]
Logical = list[Pos]
