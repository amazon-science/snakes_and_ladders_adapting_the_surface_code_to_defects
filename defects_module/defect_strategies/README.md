# Defects strategies
Module that generate the best strategies for defective rotated surface code patches. Here a
strategy is composed of both a set of superstabilizers and a valid boundary deformation.
Valid sets of superstabilizers (different possibilities) are computed first in an extended window
(a window that is bigger than the patch size such that no boundary deformation is ever needed).
After defining defect clusters (a cluster of defects that is inseparable), each of them is passed
to the module `snakes_and_ladders`, which enables the user to use both the efficient strategy and
the Auger method, to compute combinations of superstabilizers. Then, all possible and valid boundary
deformations are computed for each set of superstabilizers back in the initial window. They are
calculated with the module `boundary_deformation` which defines the holes and the different ways
to go around and cut through holes that are *truncated* back in the initial window. The effective
code distance for each strategy is computed with the module `code_distance.py`.

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)
- [Terminology](#terminology)
- [Examples](#examples)
- [Help](#help)

## Usage
The module exports different functions. The most important one is [`get_valid_strategies_in_patch`](#get_valid_strategies_in_patch)
found in `get_valid_strategies.py`. This function returns all possible and valid strategies for a
defective rotated surface code patch. The strategies are stored in [`CombinedCluster`](#CombinedCluster) data classes,
defined in `defects_clusters.py`. [*IMPORTANT: This function only works for a defective rotated surface code patch!*]
```
def get_valid_strategies_in_patch(
    patch: RotatedSurfaceCode,
    ancilla_defects: set[Pos],
    data_defects: set[Pos],
    link_defects: set[tuple[Pos, Pos]],
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
    add_padding: bool = True,
) -> list[CombinedCluster]:
    """Function that returns all the possible strategies for the defective surface code patch.
    These strategies include different choices of superstabilizers and different possible boundary
    deformations, which were all checked to be valid.

    The function works as follows:
    1. We define the *extended* and *initial* windows.
    2. We make sure that the defects configuration is valid. It happens at the second step
    because the initial window includes padding ancillas.
    3. We separate code ancilla defects from padding ancilla defects.
    4. We generate the defect clusters in the *extended* window.
    5. We run the Snakes and Ladders algorithm (we use the efficient strategy if repurpose_ancillas = True, else we use Auger's)
    in the *extended* window using the defect clusters found in 4. At this step, we use
    heuristics if any were passed. By default, runtime and memory grow exponentially
    in the number of ancilla and link defects (~ 2^n). Different heuristics can make
    this growth sub-exponential or polynomial instead.
    6. For each possible set of superstabilizers (strategy) we run the boundary deformation
    algorithm in the *initial* window: for each hole found in the *extended* window, we need
    to determine if it is closed or open in the *initial* window. For each possible boundary
    deformation, we collect all the information necessary to build a patch (i.e. a
    `CombinedCluster` object) and compute the effective code distance.
    7. We return a list of all possible strategies, i.e. a list of `CombinedCluster` objects.
    Each CombinedCluster object has fields for the effective distance, superstabilizers and
    undamaged stabilizers.

    Input arguments:
        `patch`: A RotatedSurfaceCode object which is the rotated surface code
        patch in which are found the defective components.
        `ancilla_defects`: A set of all the ancilla defects in the patch.
        `data_defects`: A set of all the data defects in the patch.
        `link_defects`: A set of all the link defects in the patch.
        `repurpose_ancillas`: An option to use snakes and ladders instead of Auger.
        The default is True.
        `use_heuristics`: Option to add heuristics when finding the different possible
        superstabilizers with snakes and ladders.
        `add_padding`: An option to add padding ancillas.

    Output arguments:
        A list of `CombinedCluster` objects corresponding each to a unique, valid strategy.
    """
```
The module exports the classes [`Heuristics`](#Heuristics) (from the `snakes_and_ladders` module and which
helps the user define the heuristics to be used for the snakes and ladders algorithm),
[`DistanceFilterFunctions`](#DistanceFilterFunctions) and [`FilterDistance`](#FilterDistance) (both defined
in `filter_strategies_by_distance.py`) which help the user sort and filter the different strategies based on
their effective code distance. Here is a minimal example:
```
heuristics = qt.Heuristics(
    minimize_zombie_qubits=False,
    minimize_max_distance_loss_per_cluster=False,
    minimize_sum_distance_loss_per_cluster=False,
)

filter_distance = qt.FilterDistance(filter_func=qt.DistanceFilterFunctions.return_best_strategies)
```
The module also exports the exception class `NoStrategyFoundError`, as well as other helper functions,
`stabilizers_commute` and `stabilizers_connected_components`. (Refer to `__init__.py` for the locations
of these functions.)

## Algorithm
Description of the algorithm to find optimal and valid strategies for a defective rotated surface code patch.

### get_valid_strategies_in_patch
This function first calls the class method [`Window.make_windows`](#Window.make_windows) defined in
`make_windows.py` in the `Window` class, to build the *initial* and *extended* windows for the patch.
The function then calls the subroutine [`generate_defects_clusters`](#generate_defects_clusters) defined
in `defects_cluster.py` to clusterize the defects. For each defects cluster, it then runs
[`run_snakes_and_ladders`](#run_snakes_and_ladders) defined in the submodule `snakes_and_ladders`.
For each strategy (set of superstabilizers) returned by snakes and ladders, then calls
[`run_boundary_deformation`](#run_boundary_deformation) defined in the submodule `boundary_deformation`.
For each strategy (set of superstabilizers AND valid boundary deformation) it then computes the effective
code distance with the function [`compute_effective_distance`](#compute_effective_distance) defined in
`code_distance.py`. Each strategy is returned as a [`CombinedCluster`](#CombinedCluster) object.

### Window.make_windows
`Window` is a class used to store data and ancilla qubits of a patch within some rectangular region
defined by a distance (the code distance of the patch as if it does not have defects) and
an offset. It also differentiates between ancillas in the code and part of the padding used
for the `efficient` strategy. The padding consists of superfluous ancilla qubits surrounding
the rotated surface code patch. They are a priori inactivem and only ever used by the
`efficient` strategy to compensate for an disabled ancilla in the patch. The `efficient`
strategy repurpose a pair of ancillas sandwiching an ancilla defect every other round to
replace the damaged stabilizer. In the case of link defects, only one neighbouring ancilla
needs to be repurposed every other round. The ckass method `make_windows` returns the initial
and extended windows for a given patch. The *extended* window is used in `snakes and ladders`
algorithm whereas the *initial* window is used in the `boundary deformation` algorithm.

### generate_defects_clusters
This function generates `DefectsCluster` objects, one for each cluster of defects in the patch.
Here a DefectsCluster object is a subpatch which encloses all defective qubits that cannot be
handled by separate super-stabilizers with the Auger strategy. This is worst case scenario that
sets an upper bound on the size of the super-stabilizeres. With the efficient strategy, we most
likely end up with smaller superstabilizers for the same cluster of defects. To build the clusters,
we start with a graph that only contains the defective qubits (here treat the qubits on both ends of
each defective link as defects). We apply the Auger method to the defects to identify the surrounding
qubits that need to be disabled or used in gauge checks, and add those qubits to the graph. The edges
in the graph represent links between qubits. After finalizing the graph, we divide it into connected
components so that each contains information of a cluster. This function is found in `defects_clusters.py`.

### run_snakes_and_ladders
A function that searches and returns a set of valid strategies for a given defect cluster. If a
heuristic is given, the function will reduce the serach space guided by the heuristic. This function
is found in the submodule `snakes_and_ladders`. Refer to the README.md there.

### run_boundary_deformation
A function that runs the whole algorithm for the boundary deformation. The function works as follows:
1. We determine which holes are open and which are closed in the initial window. We then clean the gauges
around the open holes.
2. We then get the stabilizers given the current strategy (i.e. superstabilizers defined in the extended
window) for each possible boundary deformation: we have open holes might have more than one valid deformation
due to corner placement.
This function is found in the submodule `boundary_deformation`. Refer to the README.md there.

### compute_effective_distance
A function that computes the effective distance of the deformed patch from detector graphs, given lists
of the efficient checks used and superstabilizers. The effective distance (which is the length of the
shortest logical in each direction) is based on the detector graphs with the detectors of same type
connected to a pair of virtual boundary nodes on opposite sides in the direction of the logical of
opposite type. The effective distance also depends on how the boundary is deformed. Here we include this
deformation by adding the data qubits on the deformed boundary to the graphs and attach them to the virtual
boundary nodes such that we can correctly include the distance drop between the virtual boundary nodes and
the detectors. This function is found in `code_distance.py`.

### CombinedCluster
A dataclass that contains all DefectClusters inside a window or patch, which are combined together.
The class stores the combined superstabilizers, the undamaged stabilizers in the patch and the effective
distance of the patch. This class basically has all the information needed about the strategy being
applied, and the properties needed to build a patch object. This class is found in `defects_clusters.py`.