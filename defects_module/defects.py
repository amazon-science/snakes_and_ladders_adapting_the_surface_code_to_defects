"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

import itertools

from numpy import zeros

from defects_module.base import (
    DistanceT,
    Logical,
    PauliT,
    Pos,
    Stabilizer,
    SuperStabilizer,
)
from defects_module.code_library import RotatedSurfaceCode
from defects_module.defect_strategies import (
    CombinedCluster,
    FilterDistance,
    Heuristics,
    NoStrategyFoundError,
    Window,
    find_best_subpatch_strategy,
    get_valid_strategies_in_patch,
    stabilizers_commute,
    stabilizers_connected_components,
)
from defects_module.translate import PhysicalCircuit

############################################################
# DEFECTIVE SURFACE CODE PATCH

# Patch that supports defects handled with Snakes & Ladders
############################################################


def get_frozen_qubits(patch: RotatedSurfaceCode, stabilizers: set[Stabilizer]) -> set[Pos]:
    """A data qubit in the patch is 'frozen' if it is contained in at least one
    stabilizer, and the stabilizers which contain it are only a single type.
    Returns the set of frozen data qubits in the patch.
    """
    check_matrix = zeros((len(patch.data_qubits), 2))
    for stab in stabilizers:
        for qubit in stab.data_qubits:
            check_matrix[patch.data_qubits.index(qubit), int(stab.type == PauliT.Z)] += 1
    return set(
        [
            qubit
            for i, qubit in enumerate(patch.data_qubits)
            if 0 < sum(check_matrix[i, :]) and check_matrix[i, 0] * check_matrix[i, 1] == 0
        ]
    )


def check_stabilizers_commute(stabilizers: set[Stabilizer]) -> bool:
    """Check that all stabilizers in the patch commute."""
    x_stabilizers = [stab for stab in stabilizers if stab.type == PauliT.X]
    z_stabilizers = [stab for stab in stabilizers if stab.type == PauliT.Z]
    for s1, s2 in itertools.product(x_stabilizers, z_stabilizers):
        if not stabilizers_commute(s1, s2):
            return False
    return True


def num_encoded_qubits(data_qubits: set[Pos], stabilizers: set[Stabilizer]) -> int:
    """Computes the dimension of the gauge space as
        #(data qubits) - #(stabilizers).
    Then counts the number of logical operators of the gauge space as
        #(gauge checks) - #(superstabilizers).
    Returns the first number minus half the second. Should equal 1.
    """
    gauge_dimension = len(data_qubits) - len(stabilizers)
    num_logical_operators = sum(
        [len(s.gauges) - 1 for s in stabilizers if isinstance(s, SuperStabilizer)]
    )
    if num_logical_operators % 2 == 1:
        # this is not a valid configuration
        return -1
    return gauge_dimension - num_logical_operators // 2


class DefectiveSurfaceCode(RotatedSurfaceCode):
    """Standard rotated surface code with defective qubits/links allowed.
    Defective qubits/links may not be used, and we account for this by introducing
    superstabilizers formed as a product of gauge checks.

    The choice of how to handle defects and construct superstabilizers is
    performed prior to this object, and the superstabilizers are merely given as
    inputs. It's assumed that checks are done using the same schedule as for
    the RotatedSurfaceCode, though this can almost certainly be improved.
    """

    # XXX: We have added several of the field inherited from `strategy` as separate
    #      fields to initialization to allow the user to initialize a
    #      DefectiveSurfaceCode object without a defects strategy.
    #      This can be removed if desired.

    # TODO Is it possible to make a more general "add_detectors" method?

    def __init__(
        self,
        patch: RotatedSurfaceCode,
        # Typically these fields are deduced from the strategy.
        effective_distance: DistanceT,
        undamaged_stabilizers: set[Stabilizer],
        x_superstabilizers: set[SuperStabilizer],
        z_superstabilizers: set[SuperStabilizer],
        ancilla_defects: set[Pos],
        data_defects: set[Pos],
        link_defects: set[tuple[Pos, Pos]],
        first_super_type: PauliT = PauliT.X,
        strategy: CombinedCluster | None = None,
    ) -> None:
        # initialize patch using the parent RotatedSurfaceCode `__init__` using the optimal window
        super().__init__(
            distance=patch.distance,
            sw_offset=patch.sw_offset,
            sw_type=patch.sw_type,
            vertical_logical=patch.vertical_logical,
        )

        # other patch properties related to defects
        self.ancilla_defects = ancilla_defects
        """Ancilla defects"""
        self.data_defects = data_defects
        """Data defects"""
        self.qubit_defects = self.ancilla_defects | self.data_defects
        """Dead qubits"""
        self.link_defects = link_defects
        """Link defects"""
        self.effective_distance = effective_distance
        """Effective distance of the code, e.g. (horizontal, vertical)."""
        self.window_lims = patch.extent
        """ Window limits (x0, x1, y0, y1). """
        
        # stabilizers
        self.undamaged_stabilizers = undamaged_stabilizers
        """List of undamaged stabilizers"""

        self.first_super_type = first_super_type
        """Type of superstabilizers to be checked first"""
        if first_super_type == PauliT.X:
            self.first_super = x_superstabilizers
            """First type of superstabilizers"""
            self.second_super = z_superstabilizers
            """Second type of superstabilizers"""
        elif first_super_type == PauliT.Z:
            self.first_super = z_superstabilizers
            self.second_super = x_superstabilizers

        self._stabilizers: set[Stabilizer] = (
            undamaged_stabilizers | x_superstabilizers | z_superstabilizers
        )
        """List of patch stabilizers"""

        # logicals and plotting tools
        self.logicals: dict[PauliT, Logical] = dict()
        """Logical strings along the boundaries of the patch, used for
        self.measure()."""
        self.deformed_logicals: dict[str, Logical] = dict()
        """All logical strings along the boundaries of the patch. Used in
        plotting."""
        self.efficient_checks: list = []
        """ List of efficient checks"""
        self.holes: list = []
        """List of holes"""
        if strategy:
            paulis = (
                [PauliT.X, PauliT.Z] if self.vertical_logical == PauliT.X else [PauliT.Z, PauliT.X]
            )
            for pauli, label in zip(paulis, ["left", "bottom"], strict=True):
                self.logicals[pauli] = strategy.boundaries[label]
        else:
            self.logicals = patch.logicals
        
        if strategy:
            self.deformed_logicals = strategy.boundaries
            self.detector_graphs = {
                "horizontal": strategy.horizontal_distance_graph.edges(data=True),
                "vertical": strategy.vertical_distance_graph.edges(data=True),
            }
            self.efficient_checks = strategy.efficient_strategy
            self.holes = strategy.holes

        # attributes used for constructing check circuits.
        self.ancilla_gauges = self._make_ancilla_gauges()
        """For each ancilla, hold a list of gauges checks which use it. This
        is used for constructing detectors."""

        # NOTE: The following are pieces of state which are modified as an
        # experiment is constructed. Morally, they should be contained in
        # PhysicalCircuit.
        self.check_round: int
        """The number of check rounds which have been performed after data qubit
        initialization."""
        self.gauge_msmt_inds: dict[Stabilizer, int]
        """For each gauge, track the last measurement index at which its ancilla
        was measured. This is used for constructing detectors."""

    def set_first_super_type(self, paulit: PauliT):
        """Sets the pauli type of superstabilizer that is checked immediately
        after initialization.
        """
        if paulit != self.first_super_type:
            self.first_super, self.second_super = self.second_super, self.first_super
        self.first_super_type = paulit

    def _make_ancilla_gauges(self):
        """Build a dictionary whose keys are ancillas, and whose values are lists
        of superstabilizer gauge checks containing that ancilla.
        """
        ancilla_gauges = {}
        for ss in self.first_super | self.second_super:
            for gauge in ss.gauges:
                anc = gauge.only_ancilla
                ancilla_gauges[anc] = ancilla_gauges.get(anc, []) + [gauge]
        return ancilla_gauges

    def _check(self, circ: PhysicalCircuit, stabilizers: list[Stabilizer] | None = None):
        """Perform one check round"""
        # Convert into gauges for checks
        if not stabilizers:
            stabilizers = []

        stabilizers_and_gauges = SuperStabilizer.decompose(stabilizers)
        super()._check(circ, stabilizers_and_gauges)

        # Increment the check round
        self.check_round += 1

        # For each superstabilizer we checked, shift the measurement indices of
        # of all gauges sharing ancillas with it backwards.
        for stab in stabilizers:
            if isinstance(stab, SuperStabilizer):
                for anc in stab.ancilla:
                    for gauge in self.ancilla_gauges[anc]:
                        self.gauge_msmt_inds[gauge] -= 1

    def initialize(self, circ: PhysicalCircuit, paulit: PauliT):
        """Initializes the logical patch in the given Pauli basis, then performs
        *two* rounds of checks -- one for each superstabilizer type.

        For each stabilizer that commutes with the initialization basis, we generate
        a deterministic detector.

        In the case of a defective surface code, we can be more greedy. Here, the stabilizers are
        either undamaged (the same as in the underlying surface code without defects) or
        superstabilizers.

        If the initialization basis commutes with the first type of superstabilizers we measure,
        the gauges of that superstabilizer are deterministic on the initialization plane.
        And setting up detectors using the gauges produces a lower logical error rate than
        one using the superstabilizers.
        """
        # Reset the check round and gauge measurement indices
        self.check_round = 0
        self.gauge_msmt_inds = {
            gauge: -1 for ss in self.first_super | self.second_super for gauge in ss.gauges
        }

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

        # Do first round of circuit checks
        self._check(circ=circ, stabilizers=list(self.first_super | self.undamaged_stabilizers))

        # Add detectors at undamaged stabilizers of the same Pauli type as the
        # initialization basis
        for stab in self.undamaged_stabilizers:
            if stab.type != paulit:
                continue
            circ.add_detector([(stab.only_ancilla, -1)])

        # Add detectors at gauges of the same Pauli type as the initialization basis
        for stab in self.first_super:
            if stab.type == paulit:
                for anc in stab.ancilla:
                    circ.add_detector([(anc, -1)])
            for gauge in stab.gauges:
                self.gauge_msmt_inds[gauge] = -1

        # Do second round of circuit checks
        self._check(circ=circ, stabilizers=list(self.second_super | self.undamaged_stabilizers))

        # Add detectors at every undamaged stabilizer
        for stab in self.undamaged_stabilizers:
            (anc,) = stab.ancilla
            circ.add_detector([(anc, -1), (anc, -2)])

        # For superstabilizers of the same type as the initialization basis,
        # add a detector for the whole superstabilizer
        for stab in self.second_super:
            if stab.type == paulit:
                circ.add_detector([(anc, -1) for anc in stab.ancilla])
            for gauge in stab.gauges:
                self.gauge_msmt_inds[gauge] = -1

    def measure(self, circ: PhysicalCircuit, paulit: PauliT):
        """Measure the patch in the the logical `basis`, handling superstabilizer decomposition.

        For each stabilizer that commutes with the measurement basis, we generate
        a deterministic detector.

        In the case of a defective surface code, we can be more greedy. Here, the stabilizers are
        either undamaged (the same as in the underlying surface code without defects) or
        superstabilizers.

        If the measurement basis commutes with the last type of superstabilizers we measure,
        the gauges of those superstabilizers are deterministic on the measurement plane.
        And setting up detectors using the gauges produces a lower logical error rate than
        one using the superstabilizers.
        """
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

        # Add detectors at undamaged stabilizers of the measured type
        for stab in self.undamaged_stabilizers:
            if stab.type != paulit:
                continue
            circ.add_detector([(stab.only_ancilla, -1)] + [(q, -1) for q in stab.data_qubits])

        # Collect superstabilizers of the measured type.
        first_super_is_deterministic = self.first_super_type == paulit
        deterministic_supers = (
            self.first_super if first_super_is_deterministic else self.second_super
        )

        second_super_last_round = self.check_round % 2 == 0

        if second_super_last_round ^ first_super_is_deterministic:
            # If they were *just* measured, add a detector at each gauge
            for stab in deterministic_supers:
                for gauge in stab.gauges:
                    circ.add_detector(
                        [(gauge.only_ancilla, -1)] + [(q, -1) for q in gauge.data_qubits]
                    )
        else:
            # Otherwise, add a detector at each superstabilizer
            for stab in deterministic_supers:
                circ.add_detector(
                    [(gauge.only_ancilla, self.gauge_msmt_inds[gauge]) for gauge in stab.gauges]
                    + [(q, -1) for q in stab.data_qubits]
                )

        # Add observable.
        circ.add_observable(self.logicals[paulit])

    def check(self, circ: PhysicalCircuit):
        """Performs one check round."""
        # Decide which stabilizers we'll check.
        if self.check_round % 2 == 0:
            supers = self.first_super
        elif self.check_round % 2 == 1:
            supers = self.second_super
        else:
            raise NotImplementedError

        # Do the circuit checks
        self._check(circ=circ, stabilizers=list(self.undamaged_stabilizers | supers))

        # Add a detector at each stabilizer just measured.
        for stab in self.undamaged_stabilizers:
            anc = stab.only_ancilla
            circ.add_detector([(anc, -1), (anc, -2)])

        for stab in supers:
            circ.add_detector(
                [
                    (gauge.only_ancilla, ind)
                    for gauge in stab.gauges
                    for ind in (-1, self.gauge_msmt_inds[gauge])
                ]
            )
            for gauge in stab.gauges:
                self.gauge_msmt_inds[gauge] = -1

    @classmethod
    def make_defective_patches(
        cls,
        patch: RotatedSurfaceCode,
        ancilla_defects: set[Pos],
        data_defects: set[Pos],
        link_defects: set[tuple[Pos, Pos]],
        repurpose_ancillas: bool = True,
        use_heuristics: Heuristics = Heuristics(),
        add_padding: bool = True,
        filter_distance: FilterDistance = FilterDistance(),
        first_super_type: PauliT = PauliT.X,
    ) -> tuple[int, list[DefectiveSurfaceCode]]:
        """Args:
            patch: A `RotatedSurfaceCode` patch corresponding to the initial
            patch
            ancilla_defects: A set of `Pos` corresponding to defective ancillas
            data_defects: A set of `Pos` corresponding to defective data qubits
            link_defects: A set of tuple of `Pos` corresponding to defective
            two-qubit gates in the surface code
            repurpose_ancillas: Whether we use snake n ladders or only Auger.
            use_heuristics: If we consider all orientations for the repurposing
            or only the efficient repurposing. Note that we also map link defects
            up the computation for larger code distances.
            add_paddng: A bool which allows to add padding ancillas to the patch
            The optimization would be over the effective distance. For example,
            The default value would first it maximizes the smallest distance
            along both directions, then maximizes the sum distance and returns
            the best strategies.
            first_super_type: the type of the superstabilizers that will be
            measured first.

        Returns:
            A list of `DefectiveSurfaceCode` patches for each strategy.
        """
        # compute all optimal strategies for the damaged patch
        # compute all valid strategies if any
        strategies = get_valid_strategies_in_patch(
            patch,
            ancilla_defects,
            data_defects,
            link_defects,
            repurpose_ancillas=repurpose_ancillas,
            use_heuristics=use_heuristics,
            add_padding=add_padding,
        )

        # number of strategies found
        num_strategies = len(strategies)

        if num_strategies == 0:
            raise NoStrategyFoundError("No valid strategies found for these defects.")

        # filter strategies based on their effective distance
        strategies = filter_distance.filter(strategies)

        return num_strategies, [
            cls(
                patch,
                strategy.effective_distance,
                set(strategy.undamaged_stabilizers),
                set(strategy.x_superstabilizers),
                set(strategy.z_superstabilizers),
                ancilla_defects,
                data_defects,
                link_defects,
                first_super_type=first_super_type,
                strategy=strategy,
            )
            for strategy in strategies
        ]

    def stabilizers_commute(self) -> bool:
        """Returns False if there is a pair of (super)stabilizers of the code
        which don't pairwise commute.
        """
        return all(
            [stabilizers_commute(s1, s2) for s1, s2 in itertools.combinations(self.stabilizers, 2)]
        )

    def has_frozen_qubits(self) -> bool:
        """Returns True if the patch has frozen data qubits, i.e. those acted
        upon by only stabilizers of a single Pauli type.
        """
        return bool(get_frozen_qubits(self, set(self.stabilizers)))

    def num_encoded_qubits(self) -> int:
        """Returns the number of encoded logical qubits."""
        return num_encoded_qubits(self.data_qubits, set(self.stabilizers))

    def is_connected(self) -> bool:
        """Returns False if there are two data qubits in the code which are not
        connected by a string of stabilizers.
        """
        return len(stabilizers_connected_components(set(self.stabilizers))) < 2

    def make_subpatch(self, distance: DistanceT, sw_offset: Pos):
        """Returns a DefectiveSurfaceCode patch for a region of the patch. This
        function can be used in the context of lattice surgery to define optimal
        patches to merge/split.
        """
        # make the window for the subpatch
        sw_type = self._get_stabilizer_type(sw_offset.x + 1, sw_offset.y + 1)
        undamaged_subpatch = RotatedSurfaceCode(
            distance=distance, sw_offset=sw_offset, sw_type=sw_type
        )
        window = Window.make_window(undamaged_subpatch, add_padding=True)
        window.padding_ancillas &= set(self.ancilla_qubits)

        # strategy of the subpatch
        strategy = find_best_subpatch_strategy(window, self.efficient_checks, self.holes)
        if strategy is None:
            raise NoStrategyFoundError("Could not find a valid strategy for the subpatches.")

        # defects in the the subpatch
        subpatch_qubits = window.data_qubits | window.ancillas | window.padding_ancillas
        data_defects = {q for q in self.data_defects if q in subpatch_qubits}
        ancilla_defects = {q for q in self.ancilla_defects if q in subpatch_qubits}
        link_defects = {link for link in self.link_defects if set(link).issubset(subpatch_qubits)}

        return DefectiveSurfaceCode(
                undamaged_subpatch,
                strategy.effective_distance,
                set(strategy.undamaged_stabilizers),
                set(strategy.x_superstabilizers),
                set(strategy.z_superstabilizers),
                ancilla_defects,
                data_defects,
                link_defects,
                first_super_type=self.first_super_type,
                strategy=strategy,
        )
