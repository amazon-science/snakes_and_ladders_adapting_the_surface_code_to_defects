# Make Boundaries
Module to make boundaries for a defective rotated surface code patch.
This is a submodule of `qec_toolkit.defect_strategies.boundary_deformation` and is used in the context of
boundary deformation for a defective rotated surface code patch.

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)

## Usage
The single function returned by this module is [`contour_patch`](#contour_patch) found in `contour_patch.py`.
```
def contour_patch(
    window: Window,
    boundary_deformation: BoundaryDeformation,
    undamaged_stabilizers: set[Stabilizer],
    superstabilizers: set[SuperStabilizer],
) -> tuple[set[Stabilizer], set[SuperStabilizer], dict[str, list[Pos]]]:
    """Function that returns the boundaries of the patch given the
    boundary defomration and the stabilizers. The stabilizers could be updated
    when trying to draw the boundaries: some stabilizers that fall outside
    the boundaries are removed.

    Input arguments:
        `window`: A `Window` object containing information about the initial window.
        `boundary_deformation`: A `BoundaryDeformation` object containing information
        about all the boundary deformations around the open holes in the patch.
        `undamaged_stabilizers`:  A set of all the undamaged stabilizers in the patch,
        defined in the initial window.
        `superstabilizers`: A set of all the superstabilizers in the patch, also defined
        in the initial window.

    Output arguments:
        A tuple containing:
        1) The updated set of undamaged stabilizers in the patch.
        2) The updated set of superstabilizers in the patch.
        3) A dictionary with keys as strings (`left`, `right`, `bottom` and `top`),
        and with values as lists of Pos objects, containg information about the boundary along
        each of the four sides of the patch. The list of Pos is the list of data qubits used
        to define the boundaries.
    """
```
This module also returns an exception class `NoBoundariesFoundError` found in `draw_boundaries.py`.

## Algorithm
Description of the algorithm used to define valid boundaries for the damaged patch.

### contour_patch
This is the main function of the module, found in `contour_patch.py`. It calls two subroutines:
1. [`make_boundaries`](#make_boundaries) found in `draw_boundaries.py`;
2. [`remove_extrusions`](#remove_extrusions) found in `remove_extrusions.py`.

### make_boundaries
This function creates the boundaries of the patch by returning a dictionary with the list of qubits
for each boundary (left, bottom, right, top):
1. It creates a graph for each boundary. Each graph is built from the connectivity of the stabilizers
    and gauges in the superstabilizers: each edge connects two data qubits if they are connected to the
    same ancilla qubit. The graph nodes are only data qubits that are either part of the initial
    non-deformed boundaries OR part of the boundary deformations obtained from the open holes.
2. It extracts boundaries from the graphs by first identifying possible corners, and second find the
    shortest path between each pair of those corners. We keep the corners that give the shortest paths
    overall.
3. It validates the boundaries. We make sure that:
    - No opposite boundaries must be touching.
    - Every pair of adjescent boundaries share a single data qubit.
    - There are four corners.
    - The number of X and Z checks on each data qubit is correct.
    If no valid boundaries are found, the function raises a `NoBoundariesFoundError`, also defined in
    `draw_boundaries.py`.

### remove_extrusions
This function removes all extrusions from the patch, i.e. all the stabilizers that fall outside the
boundaries that were drawn by [`make_boundaries`](#make_boundaries):
1. It defines a contour from the boundaries which is used to instantiate a `shapely.geometry.polygon.Polygon`
object.
2. It identifies all qubits that are not part of the polygon. Those are the extrusions.
3. It removes or modifies all the stabilizers that include those qubits. This means that the stabilizers of
the patch are *updated* by this function.