"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from defects_module.base import (
    PauliT,
    Pos,
)
from defects_module.code_library import (
    RotatedSurfaceCode,
)

from .utility import split_ancillas_by_native_type


class Window:
    """Class used to store data and ancilla qubits of a patch within some rectangular region
    defined by a distance (the code distance of the patch as if it does not have defects) and
    an offset. It also differentiates between ancillas in the code and part of the padding used
    for the `efficient` strategy. The padding consists of superfluous ancilla qubits surrounding
    the rotated surface code patch. They are a priori inactivem and only ever used by the
    `efficient` strategy to compensate for an disabled ancilla in the patch. The `efficient`
    strategy repurpose a pair of ancillas sandwiching an ancilla defect every other round to
    replace the damaged stabilizer. In the case of link defects, only one neighbouring ancilla
    needs to be repurposed every other round.
    """

    def __init__(
        self,
        patch: RotatedSurfaceCode,
        distance: tuple[int, int],
        sw_offset: Pos,
        add_padding: bool = True,
        sw_type: PauliT = PauliT.Z,
    ) -> None:
        """Constructor. The region of the patch used to define this object is determined by
        the distance, a tuple, and the sw_offset which is here a tuple of ints and not a Pos.
        The window encloses a reduced version of the patch with some padding (ancillas).

        Input arguments:
            `patch`: A RotatedSurfaceCode patch used to define the qubits inside the window,
            the Pauli type of vertical logical and the Pauli type of the stabilizer at the
            origin.
            `distance`: A tuple corresponding to the code distance of the patch living inside
            the window.
            `sw_offset`: The southwestern offset corner of the window, here defined as a Pos
            just like in RotatedSurfaceCode.
            `add_padding`: Option to add padding ancillas.
            `sw_type`: Pauli type of the southwestern stabilizer in the patch.
        """
        # Define the x and y limits of the region of the window
        x, y = (sw_offset.x - 1, sw_offset.y - 1)
        xlims = (x, x + 2 * distance[0])
        ylims = (y, y + 2 * distance[1])

        # store geometric properties of the window
        self.distance = distance
        """Window distance"""
        self.sw_offset = sw_offset
        """Window south-western offset"""
        self.xlims = xlims
        """horizontal limits of the window (inclusive)"""
        self.ylims = ylims
        """vertical limits of the window (inclusive)"""
        self.data_qubits = set([data for data in patch.data_qubits if self.in_window(data)])
        """List of data qubits in the window"""
        # add padding ancillas
        self.padding_ancillas: set[Pos] = set()
        if add_padding:
            # add padding on top and bottom boundaries
            self.padding_ancillas |= set(
                Pos(k, i)
                for k in range(xlims[0] + 2, xlims[1], 2)
                for i in ylims
                if Pos(k, i) not in patch.ancilla_qubits
            )
            # add padding on left and right boundaries
            self.padding_ancillas |= set(
                Pos(i, k)
                for k in range(ylims[0] + 2, ylims[1], 2)
                for i in xlims
                if Pos(i, k) not in patch.ancilla_qubits
            )
        x_ancillas, z_ancillas = split_ancillas_by_native_type(patch)
        self.x_ancillas = set(ancilla for ancilla in x_ancillas if self.in_window(ancilla))
        """List of x ancilla qubits in the window"""
        self.z_ancillas = set(ancilla for ancilla in z_ancillas if self.in_window(ancilla))
        """List of z ancilla qubits in the window"""
        self.ancillas = self.x_ancillas | self.z_ancillas
        """List of all ancilla qubits in the window"""
        self.vertical_logical = patch.vertical_logical
        """Logical on the vertical (same as the patch)"""
        self.sw_type = sw_type
        """South western ancilla type"""

    def in_window(self, qubit: Pos) -> bool:
        """Method that checks if a qubit is inside the window.

        Input arguments:
            `qubit`: A Pos corresponding to the qubit that is
            checked.

        Output arguments:
            True or False.
        """
        return (self.xlims[0] <= qubit.x <= self.xlims[1]) and (
            self.ylims[0] <= qubit.y <= self.ylims[1]
        )

    def on_boundary(self, qubit: Pos) -> bool:
        """Method that checks if the data qubit is on the boundary.

        Input arguments:
            `qubit`: A Pos corresponding to the qubit that is
            checked.

        Output arguments:
            True or False.
        """
        return (
            qubit.x == self.xlims[0] + 1  # left boundary
            or qubit.x == self.xlims[1] - 1  # right boundary
            or qubit.y == self.ylims[0] + 1  # bottom boundary
            or qubit.y == self.ylims[1] - 1  # top boundary
        )

    @classmethod
    def make_window(cls, patch: RotatedSurfaceCode, add_padding: bool = True) -> Window:
        """Class method that returns the window for a given patch.

        Input arguments:
            `patch`: A RotatedSurfaceCode object which is the rotated surface code
            patch in which are found the defective components.
            `add_padding: An option to add padding ancillas.

        Output arguments:
            A Window object.
        """
        # window for the patch
        return cls(
            patch,
            patch.distance,
            patch.sw_offset,
            add_padding=add_padding,
            sw_type=patch.sw_type,
        )

    @classmethod
    def make_windows(
        cls, patch: RotatedSurfaceCode, add_padding: bool = True
    ) -> tuple[Window, Window]:
        """Class method that returns the initial and extended windows for a given
        patch. The *extended* window is used in `snakes and ladders` algorithm
        whereas the *initial* window is used in the `boundary deformation` algorithm.

        Input arguments:
            `patch`: A RotatedSurfaceCode object which is the rotated surface code
            patch in which are found the defective components.
            `add_padding: An option to add padding ancillas.

        Output arguments:
            A tuple of two Window objects: one for the *extended* window and one for
            the *initial* window.
        """
        # initial window for boundary deformation
        initial_window = cls.make_window(patch, add_padding=add_padding)
        # extended window for Snakes n Ladders
        extended_patch = RotatedSurfaceCode(
            distance=(patch.distance[0] + 4, patch.distance[1] + 4),
            sw_offset=Pos(patch.sw_offset.x - 4, patch.sw_offset.y - 4),
            sw_type=patch.sw_type,
            vertical_logical=patch.vertical_logical,
        )
        extended_window = cls.make_window(extended_patch, add_padding=True)
        return initial_window, extended_window
