"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

from __future__ import annotations

from dataclasses import dataclass

from defects_module.base import (
    Pos,
    Stabilizer,
    SuperStabilizer,
)


@dataclass
class BoundaryDeformationString:
    """Dataclass that stores information about a boundary deformation around a single
    open hole.

    Input arguments:
        `data_qubits`: Set of data qubits in the deformed boundary around the open
        hole.
        `position`: Position of the open hole in the patch.
    """

    data_qubits: set[Pos]
    """Set of data qubits that are part of the deformed boundary string"""
    position: int
    """Position of the string along the boundary defined as an integer, i.e.
        0: Left
        1: Bottom
        2: Right
        3: Top
    """


@dataclass
class BoundaryDeformationStrings:
    """Dataclass that contains a collection of BoundaryDeformationString objects
    in the patch. Each of these objects contains information about the boundary
    deformation around each open hole in the patch. This class contains a single
    class method `merge` that takes a list of BoundaryDeformationStrings objects
    and returns a single one with all their properties merged together.

    Input arguments:
        `strings`: List of the BoundaryDeformationString objects in the patch.
    """

    strings: list[BoundaryDeformationString]
    """List of all BoundaryDeformationStrings that are in the patch."""

    @classmethod
    def merge(
        cls, boundary_strings_collection: list[BoundaryDeformationStrings]
    ) -> BoundaryDeformationStrings:
        """Method that combines all BoundaryDeformationStrings objects into a single
        BoundaryDeformationStrings object.

        Input arguments:
            `boundary_strings_collection`: List of BoundaryDeformationStrings objects
            that need to be combined into a single BoundaryDeformationStrings object.

        Output arguments:
            A BoundaryDeformationStrings object.
        """
        # create a list containing the `strings` of all the BoundaryDeformationStrings
        # objects to create and initialize a single BoundaryDeformationStrings object
        strs = [
            string
            for boundary_strings in boundary_strings_collection
            for string in boundary_strings.strings
        ]
        # return a single BoundaryDeformationStrings object
        return cls(strs)


@dataclass
class BoundaryDeformation:
    """Dataclass that stores information about all the boundary deformations around
    the open holes. The class contains a single class method `merge` that takes a list
    of `BoundaryDeformation` objects and returns a single BoundaryDeformation object
    after combining their properties.

    Input arguments:
        `superstabilizers`: Set of remaining superstabilizers in the open holes after
        cleaning the boundaries of these holes.
        `boundary_checks`: Set of stabilizers along the boundary of the holes after
        cleaning. These stabilizers correspond to a subset of the gauges of the
        superstabilizers of the holes in the extended window which became open holes
        in the initial window.
        `boundary_strings`: BoundaryDeformationStrings object associated with the
        open hole.
        `deformed_corners`: List of initial and final corners in tuple format, i.e.
        (initial corner, final corner).
    """

    superstabilizers: set[Stabilizer | SuperStabilizer]
    """Set of modified superstabilizers in the open holes after cleaning the boundaries."""
    boundary_checks: set[Stabilizer]
    """Checks along the boundaries after cleaning.

    They are a subset of the gauges of the superstabilizers of the open holes defined
    in the extended window.
    """
    boundary_strings: BoundaryDeformationStrings
    """BoundaryDeformationStrings around the open holes."""
    deformed_corners: list[tuple[Pos, Pos]]
    """List of (initial, final) corners."""

    @classmethod
    def merge(
        cls, boundary_deformation_collection: list[BoundaryDeformation]
    ) -> BoundaryDeformation:
        """Method that combines all BoundaryDeformation objects into a single
        BoundaryDeformation object.

        Input arguments:
            `boundary_deformation_collection`: List of BoundaryDeformation objects
            that need to be combined into a single BoundaryDeformation object.

        Output arguments:
            A BoundaryDeformation object.
        """
        # create a set containing superstabilizers of all the BoundaryDeformation
        # objects to create and initialize a single BoundaryDeformation object
        superstabilizers = set(
            stab for el in boundary_deformation_collection for stab in el.superstabilizers
        )
        # create a set containing the boundary checks of all the BoundaryDeformation objects
        boundary_checks = set(
            stab for el in boundary_deformation_collection for stab in el.boundary_checks
        )
        # create a BoundaryDeformationStrings object from the boundary_strings of all
        # the BoundaryDeformation objects
        boundary_strings = BoundaryDeformationStrings.merge(
            [el.boundary_strings for el in boundary_deformation_collection]
        )
        # combine the deformed corners of all the BoundaryDeformation objects
        deformed_corners = [
            corner for el in boundary_deformation_collection for corner in el.deformed_corners
        ]
        # return a single BoundaryDeformation object
        return cls(
            superstabilizers,
            boundary_checks,
            boundary_strings,
            deformed_corners,
        )
