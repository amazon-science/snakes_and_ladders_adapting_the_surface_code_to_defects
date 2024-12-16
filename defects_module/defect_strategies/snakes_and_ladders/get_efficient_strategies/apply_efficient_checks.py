"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import (
    PauliT,
    Pos,
)

from .define_efficient_checks import EfficientRepurposing, TriangularCheck
from .efficient_defect_clusters import EfficientDefectsCluster


def _efficient_repurposing_per_defect(
    cluster: EfficientDefectsCluster,
    anc1: Pos,
    defect: Pos,
    anc2: Pos,
) -> list[EfficientRepurposing]:
    """Function that returns the best repurposing strategies. We try to build pairs
    of efficient checks in both directions around the ancilla of the damaged weight-4
    check defined here as `defect`. `anc1` and `anc2` are the two ancillas used
    to do the weight-4 check with two weight-2 triangular (efficient) checks that are
    *symmetric*, i.e. `anc1` and `anc2` are different from `defect`. This function
    attempts to use asymmetric checks for links when possible. It returns a list of
    all possible pairs of efficient checks (i.e. triangular data--ancilla--data checks).

    The function works as follows:
        1) We find the data qubits that are in the triangular checks.
        2) We try to make an asymmetric check if the `defect` is an ancilla zombie
        (i.e. it is not dead). A check only works if the data qubits and the ancilla
        being used are not defective.
        3) We then try to build TriangularCheck objects for symmetric checks: we
        repurpose one or two ancillas depending on whether or not we were able to
        make a valid asymmetric check in 2).
        4) We return all possibilities of efficient checks.

    Input arguments:
        `defect`: An ancilla defect or an ancilla zombie in the case of a link defect.
        `anc1`: A repurposed ancilla used to do a weight-2 triangular (efficient) check.
        `anc2`: A repurposed ancilla used to do a weight-2 triangular (efficient) check.

    Output arguments:
        A list of EfficientRepurposing objects corresponding to the different ways we
        can repurpose ancillas to reconstruct the damaged stabilizer with efficient
        weight-2 checks whend possible. Gauge checks found in the Auger method are not
        stored in that list.
    """
    # Identify the orientation of the checks
    if anc1.x == anc2.x:
        # △ check and its mirrored check
        # make sure the ancillas are ordered correctly
        if anc2.y < anc1.y:
            anc1, anc2 = anc2, anc1
        # Data qubits in the broken stabilizer, grouped
        # in pairs corresponding to the weight-2 checks
        data_qubits = [
            (defect + Pos(1, -1), defect + Pos(-1, -1)),
            (defect + Pos(1, 1), defect + Pos(-1, 1)),
        ]
    else:
        # ◁ check and its mirrored check
        # make sure the ancillas are ordered correctly
        if anc2.x < anc1.x:
            anc1, anc2 = anc2, anc1
        # Data qubits in the broken stabilizer, grouped
        # in pairs corresponding to the weight-2 checks
        data_qubits = [
            (defect + Pos(-1, 1), defect + Pos(-1, -1)),
            (defect + Pos(1, 1), defect + Pos(1, -1)),
        ]

    def make_triangle(anc: Pos, dat: tuple[Pos, Pos]) -> TriangularCheck:
        # Build a trianagular check
        data1, data2 = dat
        return TriangularCheck(anc, defect, data1, data2)

    def check_triangle(check: TriangularCheck) -> bool:
        """Check if a triangular check is feasible in this cluster."""
        # Check if one or more of its edges are defective
        dead_links = any(
            [set([edge, edge[::-1]]) & cluster.links_defective for edge in check.edges]
        )
        # Check if any qubit in the check is defective
        dead_qubits = not set(check.nodes).isdisjoint(cluster.defects)
        # Return False if the check is dead otherwise True
        if dead_links or dead_qubits:
            return False
        return True

    # These are the bottom or left check ancillas
    first_triangles: list[TriangularCheck] = []
    # These are the top or right check ancillas
    second_triangles: list[TriangularCheck] = []
    # we look at the inward triangular checks
    if defect not in cluster.ancillas_defective:
        # inward triangle for the bottom or left check
        check = make_triangle(defect, data_qubits[0])
        if check_triangle(check):
            first_triangles.append(check)
        # inward triangle for the top or right check
        check = make_triangle(defect, data_qubits[1])
        if check_triangle(check):
            second_triangles.append(check)
    # we look at the bottom or left outward triangular check
    # we only do it for link defects if we couldn't find
    # a valid asymmetric check. We always do it for ancilla
    # defects.
    if (
        not first_triangles
        and anc1 not in cluster.ancillas_defective | cluster.unavailable_padding_ancillas
    ):
        check = make_triangle(anc1, data_qubits[0])
        if check_triangle(check):
            first_triangles.append(check)
    # we look at the top or right outward triangular check
    if (
        not second_triangles
        and anc2 not in cluster.ancillas_defective | cluster.unavailable_padding_ancillas
    ):
        check = make_triangle(anc2, data_qubits[1])
        if check_triangle(check):
            second_triangles.append(check)

    # If there is no valid repurposing, return nothing
    if not first_triangles and not second_triangles:
        return [
            EfficientRepurposing(defect, []),
        ]
    # If there is a valid repurposing on only one side
    # then set the other side to have the standard repurposing
    # (even if it's not going to work on that side)
    if not first_triangles:
        return [EfficientRepurposing(defect, [check]) for check in second_triangles]
    if not second_triangles:
        return [EfficientRepurposing(defect, [check]) for check in first_triangles]

    # We combine the first and second checks to return valid repurposing strategies
    return [
        EfficientRepurposing(defect, [check1, check2])
        for check1 in first_triangles
        for check2 in second_triangles
        if check1.ancilla != check2.ancilla
    ]


def efficient_repurposing_per_cluster(
    cluster: EfficientDefectsCluster,
) -> list[list[EfficientRepurposing]]:
    """Function that pairs repurposed ancillas for the efficient strategy
    around each bulk ancilla defect in both horizontal or vertical orientations.
    For link defects, we also allow for asymmetric repurposing while
    reusing the ancilla in the defective link (unless it is defective).

    The function works as follows:
        1) For each ancilla defect / zombie in the cluster, we find all
        valid repurposings that vertical (△ and ▽) around the defect / zombie
        and all those that are horizontal (◁ and ▷).
        2) We put vertical / horizontal first if vertical / horizontal minimize
        the distance loss locally for the defect / zombie. We add them
        the list of efficient checks.

    Input arguments:
        `cluster`: A EfficientDefectsCluster object which stores information
        about the strategy such as the efficient checks being used.

    Output arguments:
        A list of all possible repurposing (i.e. a list of EfficientRepurposing objects),
        for each ancilla defect or zombie.
    """
    # Define the list of repurposed weight-2 checks
    efficient_checks: list[list[EfficientRepurposing]] = []
    # Determine which stabilizers need to be repaired
    replaced_ancillas = cluster.ancillas_defective | cluster.ancillas_zombie
    # For each broken stabilizer, we determine the different possible ways
    # to repurpose neighboring ancillas
    for defect in replaced_ancillas:
        # vertical △ and ▽ checks
        ver_checks = _efficient_repurposing_per_defect(
            cluster,
            defect + Pos(0, -2),
            defect,
            defect + Pos(0, 2),
        )
        # horizontal ◁ and ▷ checks
        hon_checks = _efficient_repurposing_per_defect(
            cluster,
            defect + Pos(-2, 0),
            defect,
            defect + Pos(2, 0),
        )
        # next we merge the lists of horizontal and vertical checks
        # first we determine the number of weight-2 checks that were found in total for the
        # vertical and horizontal orientations
        n_ver_checks = len([check for repurposing in ver_checks for check in repurposing.checks])
        n_hon_checks = len([check for repurposing in hon_checks for check in repurposing.checks])
        if n_ver_checks == n_hon_checks:
            # If we found the same number, then we put the efficient checks first followed by
            # the lossy ones to preserve the distance as much as possible
            # we determine the type of the ancilla in the patch and then determine which
            # orientation is distance-preserving
            ancilla_type = PauliT.X if defect in cluster.x_ancillas else PauliT.Z
            efficient_is_vertical = ancilla_type == cluster.vertical_logical
            checks = ver_checks + hon_checks if efficient_is_vertical else hon_checks + ver_checks
        else:
            # Otherwise we will put the orientation with the hightest number of weight-2 checks
            # first to minimize the number of zombie qubits.
            checks = (
                ver_checks + hon_checks if n_ver_checks > n_hon_checks else hon_checks + ver_checks
            )
        efficient_checks.append(checks)
    return efficient_checks
