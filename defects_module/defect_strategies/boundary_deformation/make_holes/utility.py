"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from defects_module.base import Pos


def any_on_boundary(
    lims: tuple[int, int],
    coord: str,
    qubits: set[Pos] | list[Pos],
) -> bool:
    """Helper function that checks if the contour lies on the boundary of the region
    defined by the tuple `lims` along the coordinate `coord` (x or y for horizontal
    or vertical).

    Input arguments:
        `lims`: Limits of the window along the coordinate `coord`.
        `coord`: Coordinate (x or y) of the qubits.
        `qubits`: Qubits to check.

    Output arguments:
        True or False.
    """
    return any([getattr(qubit, coord) in (lims[0] + 1, lims[1] - 1) for qubit in qubits])


class InvalidHoleFoundError(Exception):
    """Error that is raised when a data qubit or an ancilla qubit is on the boundary
    of a patch or a window. The code is then considered to be uncorrectable.
    """

    pass
