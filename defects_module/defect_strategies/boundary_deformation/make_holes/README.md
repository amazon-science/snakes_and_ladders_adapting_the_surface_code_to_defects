# Make holes
Module to make holes in a defective rotated surface code patch.
This is a submodule of `qec_toolkit.defect_strategies.boundary_deformation` and is used in the context of
boundary deformation for a defective rotated surface code patch.

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)

## Usage
The main class in this module is [`Hole`](#Hole) defined in `make_holes.py`.
```
class Hole:
    def __init__(
        self,
        window: Window,
        x_superstabilizers: list[SuperStabilizer],
        z_superstabilizers: list[SuperStabilizer],
        data_defects: set[Pos],
        ancilla_defects: set[Pos],
        ancilla_zombies: set[Pos],
    ) -> None:
        """Class that stores information about a given hole in the defective patch.
        Each hole corresponds to a pair of superstabilizers that are defined such as
        to `patch` the defective components. A hole that is closed in the initial window
        (i.e. can be patched with Auger or efficient strategies) will be kept as it is in
        the patch. However, an open hole, due to the defects being on or close to the
        boundary, will deform the patch boundary. The checks that are kept along the
        deformed boundary are computed here.

        The constructor works as follows:
        1) We compute the modified superstabilizers in the *initial window*.
        2) We check if the hole's superstabilizers are still `connected` (i.e touching)
        in the *initial window*.
        3) We check if the modified superstabilizers still commute or not.
        4) The hole is said open if the superstabilizers are no longer `connected` or
        if they don't commute anymore.
        5) We find all the `contours` (closed cycles in a graph whose connectivity is
        determined by the gates in the superstabilizers) in the hole. They help us
        compute the possible boundary deformations.
        6) We determine which defects are contained inside the hole.
        7) If the hole is open then we break the superstabilizers into smaller
        superstabilizers and stabilizers (previously gauge checks) and keep only the
        checks along the `contours` of the hole (only the contours that are found to
        be open in the *inital window*) that match the type of the boundary. We refer
        to this step as `cleaning the boundary`.

        Input arguments:
            `window`: A `Window` object containing information about the initial window.
            `x_superstabilizers`: List of X superstabilizers associated with the hole.
            `z_superstabilizers`: List of Z superstabilizers associated with the hole.
            `data_defects`: List of data defects in the patch.
            `ancilla_defects`: List of ancilla defects in the patch.
            `ancilla_zombies`: List of ancilla zombies in the patch.
        """
```
The module also contains an exception class `InvalidHoleFoundError` defined in `utility.py`.

## Algorithm
Description of the algorithm used to define holes in the defective patch and to clean
the boundaries if they are open.

### Hole
The entire algorithm takes place in the `Hole` class, found in `make_holes.py`.
Its constructor calls several subroutines. The most important ones are
1. [`pair_superstabilizers`](#pair_superstabilizers) function found in `make_holes.py`.
This function is used pair superstabilizers together to define the holes.
2. [`Cycle`](#Cycle) class found in `make_cycles.py`. Each hole has at lease one `cycle`
in it. The perimeter of the cycles are defined by the gauges of the superstabilizers of
the hole.
3. [`clean_gauges_if_open`](#clean_gauges_if_open) found in `clean_gauges_in_open_hole.py`.
This function cleans the gauges in a hole that is considered *open* in the *initial window*
(because one or many of its cycles are truncated in that window). Since the superstabilizers
of the hole no longer commute, we break them apart and remove the gauge checks of a single
Pauli type such that these open cycles are now part of the definition of the boundaries.

### pair_superstabilizers
It is a helper function that pairs superstabilizers associated with the same hole.
The function works as follows:
    1. We build a networkx graph whose nodes are the superstabilizers.
    2. For each pair of X and Z superstabilizers we check if their gauge checks
    commute or not. If any pair of gauge checks if found to not commute, then the
    pair of superstabilizers is part of the same hole and is therefored `paired`.
    3. If the pair of superstabilizers if paired, then we add an edge between
    the two superstabilizers in the graph.
    4. We find all the connected components in the graph. Each connected component
    contains all the superstabilizers that are in the same `hole`.
    We return a list of tuples (X superstabilizers, Z superstabilizers) where each
    element in the kth tuple is a list of X or Z superstabilizers contained in the
    kth hole.

### Cycle
It is a dataclass that stores information about a cycle in a hole. A `Cycle` is
defined by a cycle in a graph whose vertices and edges are determined by the
superstabilizers in the hole in the *extended window*. Every edge correspond to
a gate in one of the superstabilizers. The cycle must contain a defective qubit
(data or ancilla defect) or a defective link (ancilla zombie). The cycle is
defined by a perimeter which is a list of data qubits. The class also contains
the gauges along the cycle in the *initial window*. It also stores whether or not
the cycle is closed in the *initial window*. It also contains a `Polygon` object
whose contour is the cycle. This polygon is used to determine the qubits that
are contained or not inside the cycle.

### clean_gauges_if_open
This function returns the gauges in an open hole that are kept along the deformed
boundary that goes around or through that hole. The function works as follows:
    1. We determine if the hole is at a corner or if it's along an edge of the patch.
    2. If it is a corner hole, we call [`clean_gauges_at_corner`](#clean_gauges_on_edge).
    3. If it is a corner hole, we also check if the hole touches to the edges of the
    extended window. If it does, then we also call [`clean_gauges_on_edge`](#clean_gauges_on_edge)
    for both Pauli types: this is because it makes sense to define the new corner on both
    edges of the patch.  If it does not, then we call [`clean_gauges_on_edge`](#clean_gauges_on_edge)
    only for the Pauli type for which the new corner makes sense. Otherwise, we would
    end up with more than four corners.
    4. If it is not a corner hole but instead an edge hole, we determine the type
    of the edge on which is the hole. We then call [`clean_gauges_on_edge`](#clean_gauges_on_edge).
    Unlike the corner hole case, there is always a single BoundaryDeformation object for
    this case.
    5. We return the list of all possible boundary deformations.

#### clean_gauges_on_edge
This function is found in `clean_gauges_on_edge.py`. This function cleans gauges of the wrong
pauli type in an open hole along the edge of the patch. This function is designed to *minimize*
the amoung of digging in the patch. The function works as follows:
    1. We get the position of the hole in the patch, i.e. along which edge
    (left, bottom, right or top) it is.
    2. We find all the cycles in the hole that are 'open', i.e. truncated at
    the boundary.
    3. For each open cycle, we collect which gauges to keep and which to
    gauges to remove along that cycle. We also collect the data qubits
    along that cycle since they will be part of the deformed boundary.
    4. For each closed cycle, we also collect the gauges along that cycle.
    They will form smaller superstabilizers.
    5. We return the boundary deformation.


#### clean_gauges_at_corner
This function is found in `clean_gauges_at_corner.py`. This function cleans the gauge
checks of the wrong pauli type in an open corner hole for each possible corner position.
The function works as a follows:
    1. For each cycle within the hole, we check if the cycle touches the boundary.
    If not, we move to the next cycle.
    2. We then check that the new corner is within the initial window. We then remove
    it from the `gauge_graph` of the `Cycle` object such that we can find two connected
    components in that graph (because we removed the corner connecting them). These
    two components correspond to two boundary deformations along two perpendicular edges
    of the patch. If we found more than two connected components, we move to the next cycle.
    3. For these two connected components, we evaluate which one could be of vertical pauli
    type and which could be of horizontal pauli type. They could be interchangeable, depending
    on the hole.
    4. For each gauge checks in `gauges`, we determine if the gauge check falls in the vertical
    or horizontal boundary. We then determine if the gauge needs to be cleaned. If it is cleaned,
    then we also update the data qubits in the new boundary by adding those contained in the
    cleaned gauge check.
    5. We return the boundary deformations for each possibility (each corner, each combination of
    vertical/horizontal placement).