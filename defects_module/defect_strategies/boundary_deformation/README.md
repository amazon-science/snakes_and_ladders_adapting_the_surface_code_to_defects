# Boundary deformation
Module to generate all possible valid boundary deformations for a defective rotated surface
code patch. This is a submodule of `qec_toolkit.defect_strategies` and is used in the context of
finding an optimal set of stabilizers for a defective rotated surface code patch.

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)

## Usage
The single function returned by this module is [`run_boundary_deformation`](#run_boundary_deformation)
found in `run_boundary_deformation.py`.
```
def run_boundary_deformation(
    window: Window, cluster_combination: Sequence[DefectsCluster]
) -> tuple[list[TriangularCheck], list[CombinedCluster]]:
    """Function that runs the whole algorithm for the boundary deformation.

    The function works as follows:
        1) We determine which holes are open and which are closed in the initial window.
        We then clean the gauges around the open holes.
        2) We then get the stabilizers given the current strategy (i.e. superstabilizers
        defined in the extended window) for each possible boundary deformation: we have
        open holes might have more than one valid deformation due to corner placement.

    Input arguments:
        `window`: Window object containing information about the *initial window*.
        `cluster_combination`: Sequence of `DefectsCluster` objects from which we can
        determine the superstabilizers and the holes. These DefectsCluster objects
        are defined in the *extended* window.

    Output arguments:
        A tuple containing:
        1) A list of `TriangularCheck` objects (efficient checks) in the patch
        2) A list of `CombinedCluster` objects, one for each possible strategy
        associated with a specific choice of boundary deformation. The CombinedCluster
        object gives all the information necessary to build the patch.
    """
```

## Algorithm
Description of the algorithm used to find optimal and valid boundary deformations. The general
idea of the algorithm is as follows:
1. The superstabilizers were first obtained in the *extended* window of the damaged patch using
for example the `snakes_and_ladders` module. We use those superstabilizers to define the holes in
the patch, still in the *extended* window. We define holes using the `make_holes` submodule.
2. Using the *initial* (truncated) window we identifty which of these holes are close and which
are now open in the *initial* window. Open holes will deform the boundaries. This is determined
inside the `make_holes` submodule.
3. We deform the boundaries by first cleaning the gauges in the open holes. The idea here is to
minimize the amount of digging in the patch. In the case of a corner hole, we also need to try
different new corner placements. This is still happened inside the `make_holes` submodule.
4. Once we identify all the possible ways of deforming the boundaries *locally* around or through
the open holes, we then need to drawn the *entire* boundaries of the patch. This is now 
happening in the `make_boundaries` submodule. For that step, we need to draw valid paths connecting
the corners that maximize the size of the patch while avoiding *frozen* qubits and extrusions that
do not help the code distance but hurt the logical performance by adding more noise.

### run_boundary_deformation
This function works as follows:
    1. We create the `Hole` objects using the `Hole` class in the `make_holes` module.
    These objects hold information about the superstabilizers and all the possible
    ways to deform the boundaries around and through the open holes.
    2. We obtain all this information via a function called [`get_strategies_from_holes`](#get_strategies_from_holes)
    defind in `run_boundary_deformation.py`.

### get_strategies_from_holes
Defined in `run_boundary_deformation.py`, this function:
    1. Calls `collect_data_from_holes` (in `get_patch-stabilizers.py`) to collect all the superstabilizers,
    defects and possible *local* boundary deformations from the holes in the patch.
    2. Calls `update_patch_stabilizers` (in `get_patch-stabilizers.py`) to update the stabilizers of the patch
    by first removing all disconnected stabilizers, updating stabilizers that involve a repurposed ancilla to
    superstabilizers (this is *only* from the algorithm perspective because they are measured every other round).
    3. Calls `contour_patch` from the `make_boundaries` submodule to draw the entire boundaries of the patch.
    4. It returns a list of all possible strategies (boundary deformations) with all the properties necessary
    to build a patch for each.