# Snakes and Ladders
Module to generate all possible combinations of superstabilizers for a defective rotated surface
code patch, including efficient checks (i.e. weight-2 checks using repurposed ancillas).
This is a submodule of `qec_toolkit.defect_strategies` and is used in the context of
finding an optimal set of stabilizers for a defective rotated surface code patch.

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)

## Usage
The single function returned by this module is [`run_snakes_and_ladders`](#run_snakes_and_ladders)
found in `run_snakes_and_ladders.py`.
```
def run_snakes_and_ladders(
    defect_cluster: DefectsCluster,
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
) -> list[EfficientDefectsCluster]:
    """Function that searches and returns a set of valid superstabilizers for a given
    defect cluster. The function will reduce the search space using any heuristic that
    is passed. If repurpose_ancillas is set to True the function will use efficient
    checks (weight-2 checks using repurposed ancillas), otherwise it will run Auger.

    The function works as follows:
        1) We call the function `find_efficient_strategies_per_cluster` from the `get_efficient_strategies`
        module to find all possible strategies. Refer to the README in this module for
        more details on the algorithm.
        2) We then filter the strategies such as to keep only the unique ones. To do
        so, we compare the list of efficient checks of each strategy, and their set of
        superstabilizers.

    Input arguments:
        `defect_cluster`: A defect cluster object defined in the *extended window*.
        There is no notion of boundary deformation in this module. Any defect cluster
        is defined in a window such that holes defined by the Auger strategy will be
        completely closed. 
        `repurpose_ancillas`: Boolean which is set to True by default such that we
        allow efficient checks. If set to False, we only use the Auger strategy.
        `use_heuristics`: A Heuristics object which contains the information about
        which heuristics to use for the search of valid strategies.

    Output arguments:
        A list of EfficientDefectsCluster objects where each object contains information
        about a valid strategy (e.g. superstabilizers) for the defects cluster.
    """
```

## Algorithm
Description of the algorithm used to find valid set of superstabilizers for a given defects
configuration.

### run_snakes_and_ladders
The function works as follows:
1. We call the function `find_efficient_strategies_per_cluster` from the `get_efficient_strategies`
module to find all possible strategies. Refer to the README in this module for
more details on the algorithm.
2. We then filter the strategies such as to keep only the unique ones. To do
so, we compare the list of efficient checks of each strategy, and their set of
superstabilizers.

### set_superstabilizers

This function finds and assigns the superstabilizers of a defect cluster. The function
works as follows:
1. Define all the gauges in the defects cluster with the function [`define_checks`](#define_checks).
2. Combine the gauges into superstabilizers with the function [`combine_checks`](#combine_checks).
3. Store the superstabilizers in cluster.

### define_checks

This function lists all X-type and Z-type gauge checks in the defect cluster. The function
works as follows:
1. For every efficient check (stored as a TriangularCheck object) we define a stabilizer of
whose Pauli type correspond to the Pauli type of the damaged stabilizer that it corresponds to.
2. For every functioning ancilla (that is not however not a zombie ancilla) we define a gauge
check using all its active neighboring data qubits. The Pauli type is that of the corresponding
damaged stabilizer.

### combine_checks

This function combines the X-type and Z-type gauge checks into superstabilizers based
on their proximity to the defects. The function works as follows:
1. We define two networkx graphs: one for the X-type superstabilizers and
one for the Z-type superstabilizers.
2. For each data qubit that is either defective or a zombie data qubit we
extract all the neighbouring ancilla qubits (independently of their status)
and add an edge between the data qubit and each of these ancilla in the
graph of the Pauli type that corresponds to that of the stabilizer that
the ancilla is associated with in the undamaged patch.
3. For each efficient check (TriangularCheck), we add an edge between the
replaced ancilla and replaced ancilla in the graph of the Pauli type that
corresponds to that of the replaced ancilla in the undamaged patch. We
also add an edge perpendicular to the efficient checks, connecting the
weight-4 stabilizers in weight-8 supercheck, in the graph of the opposite
Pauli type.
4. We build a superstabilizer from each connected component in each graph.
The Pauli type of the superstabilizer is that of the graph.