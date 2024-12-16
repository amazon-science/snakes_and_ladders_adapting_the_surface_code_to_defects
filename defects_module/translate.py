"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from typing import Sequence

from stim import Circuit

from defects_module.base import GateT, Logical, NoiseModel, Pos


class PhysicalCircuit:
    """Stim circuit writer."""

    noise_model: NoiseModel
    """The noise model associated to this circuit."""
    stim_circuit: Circuit
    """The stim circuit for this experiment."""
    msmt_index: int
    """The total number of physicalmeasurements in the experiment."""
    qubit_measurement_map: dict[Pos, list[int]]
    """Mapping from qubit positions to the list of absolute measurement indices
    at which the qubit was measured (appearing in-order)."""
    observables: list[list[int]]
    """An observable is defined by a list of measurement indices. This is a list
    of observables in the experiment."""

    def __init__(self, noise_model: NoiseModel | None = None):
        if not noise_model:
            noise_model = dict()
        self.noise_model = noise_model
        self.stim_circuit = Circuit()
        self.msmt_index = 0
        self.qubit_measurement_map = {}
        self.observables = []

    def add_gate(self, gate: GateT, qubits: Sequence[Pos]):
        """Adds a gate operation with noise to the circuit."""
        noise_params = self.noise_model.get(gate, [])

        if gate == "M" and any(noise_type == "M" for noise_type, _ in noise_params):
            # TODO Fix measurement noise.
            _, noise_strength = [n for n in noise_params if n[0] == "M"][0]
            for q in qubits:
                self.stim_circuit += Circuit(f"M({noise_strength}) {int(q)}")
                self._handle_msmt_inds(q)
            return

        self.stim_circuit.append(gate, qubits)
        qubit_str = " ".join(str(int(q)) for q in qubits)
        for noise_type, noise_strength in noise_params:
            if noise_strength == 0.0:
                continue
            self.stim_circuit += Circuit(f"{noise_type}({noise_strength}) " + qubit_str)

    def _handle_msmt_inds(self, q: Pos):
        """Incrementing the measurement indices."""
        current_msmt_inds = self.qubit_measurement_map.get(q, [])
        self.qubit_measurement_map[q] = current_msmt_inds + [self.msmt_index]
        self.msmt_index += 1

    def add_detector(self, targets: list[tuple[Pos, int]]):
        """Add detector instructions, given a list of pairs (qubit, timestamp).
        For example, the pair [(q, -1), (q, -2)] means to XOR the last two
        measurement outcomes from qubit q.
        """
        det_str = " ".join(
            f"rec[{self.qubit_measurement_map[q][ind] - self.msmt_index}]" for (q, ind) in targets
        )
        self.stim_circuit += Circuit(f"DETECTOR {det_str}")

    def add_observable(self, logical: Logical):
        """Add an observable instruction, given a set of measured qubits."""
        obs_str = " ".join(
            f"rec[{self.qubit_measurement_map[q][-1] - self.msmt_index}]" for q in logical
        )
        self.stim_circuit += Circuit(f"OBSERVABLE_INCLUDE({len(self.observables)}) {obs_str}")

        self.observables.append([self.qubit_measurement_map[q][-1] for q in logical])
