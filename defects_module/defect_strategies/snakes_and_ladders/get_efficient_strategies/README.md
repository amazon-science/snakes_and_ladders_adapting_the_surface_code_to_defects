# Get Efficient Strategies
Module used to find all valid combinations of superstabilizers for a defect cluster.
This is a submodule of `qec_toolkit.defect_strategies.snakes_and_ladders`. 

## Table of Contents
- [Usage](#usage)
- [Algorithm](#algorithm)

## Usage
The single function returned by this module is [`find_efficient_strategies_per_cluster`](#find_efficient_strategies_per_cluster)
found in `find_efficient_strategies.py`: 
```
def find_efficient_strategies_per_cluster(
    cluster: EfficientDefectsCluster,
    repurpose_ancillas: bool = True,
    use_heuristics: Heuristics = Heuristics(),
) -> list[EfficientDefectsCluster]:
    """Function that iterates over all subsets of the link defects in the given cluster,
    i.e. for each link defect we disable either the neighboring ancilla qubit or the data
    qubit. For each of these configurations, we apply the Snakes & Ladders algorithm to
    find all valid strategies or reduce the search space according to the given heuristic.
    We return the strategies found for each configuration.

    The function works as follows:
        1) Create all possible assignments for the link defects to data defects or
        zombie ancillas.
        2) For each assignment, we call `_find_efficient_strategies_per_link_assignment` above
        to find all valid strategies.
        3) We filter the strategies to keep the most optimal ones if heuristics were
        passed.

    Input arguments:
        `cluster`: An EfficientDefectsCluster object that has information about
        the defect configuration.
        `repurpose_ancillas`: A boolean that determines if we use the efficient
        strategy (True) or not (False) in which case we use Auger only. The
        default value is True.
        `use_heuristics`: A Heuristics object that contains information about which
        heuristics to use.

    Output arguments:
        A list of EfficientDefectsCluster objects corresponding to the different
        strategies for the defects cluster.
    """
```
There are also a few important classes in this module:
[`Heuristics`](#Heuristics), [`TriangularCheck`](#TriangularCheck) and [`EfficientDefectsCluster`](#EfficientDefectsCluster).

## Algorithm
Description of the Snakes and Ladders algorithm used to find valid superstabilizers in the extended window of a damaged patch.

### find_efficient_strategies_per_cluster
The first step is to make the link defects to either data defects or ancilla zombies. Depending on the heuristics that is passed,
the link defects can also be restrained to ancilla zombies only to speed up the algorithm. This could however result in a smaller
effective code distance. For each assignment of link defects (it can also be an empty list if there is no link defect), the
function calls [`_find_efficient_strategies_per_link_assignment`](#_find_efficient_strategies_per_link_assignment) which will find all
the valid strategies for that assignment using heuristics if any are passed. Heuristics and information about the optimal strategies
are stored in the class [`Scores`](#Scores). The strategies are filtered to keep only the best ones if heuristics were passed.

### _find_efficient_strategies_per_link_assignment

A function that finds maximal ways of reviving zombie data qubits, where maximal is defined as using the efficient strategy over
the Auger method wherever it is possible, and generates all possible combinations of efficient checks that are valid, i.e. such that
any ancilla is only repurposed once. One of its input arguments is a boolean `use_snakes_and_ladders` which determines if we use the
efficient strategy (True) or not (False) in which case we use Auger only. The default value is True. This function knows about the
heuristics because it asks for a [`Scores`](#Scores) class as another input argument. If a heuristic is given, it narrows the search
to a subset of the combinations. If repurpose_ancillas is set to False, just return the result of Auger's algorithm. The function
works as follows:
1. If repurpose_ancillas = False, then we run the Auger method. To do so, we zombify all data qubits around the defective
or zombie ancillas. We also zombify data qubits in weight-1 checks.
2. If repurpose_ancillas = True, then we find all possible combinations of efficient checks with efficient_repurposing_per_cluster.
For each possible combination of efficient checks we zombify the data qubits that around ancilla defects/zombies that are not part of
efficient checks. We also zombify all data qubits that are part of weight-1 checks. We decide if we keep the strategy or not based on
the heuristics.

### Scores
A class that helps apply the heuristic to strategies per cluster and keeps track of optimal strategies that were found before.
It takes a [`Heuristics`](#Heuristics) object as input.

### Heuristics
```
class Heuristics:
    """Class that defines what heuristics to use to speed up the
    computation of the best strategies. The default is to have no
    heuristics at all.

    Input arguments:
        `link_defects_to_ancilla`: Whether or not we map all link
        defects to ancilla zombies. The default is False. If True,
        then link defects are never mapped to data defects. If False,
        then link defects are mapped to both data defects and ancilla
        zombies.
        `n_zombie`: The number of configurations with the same number
        of functioning data qubits being recorded before we stop the
        search. If a strategy is found to have more functioning data
        qubits then the counter is restarted. This heuristic uses the
        fact that we order the strategies such that efficient checks
        that are not lossy (that locally minimizes the distance loss)
        are tried first. This is a strong heuristic which makes the
        runtime polynomial but can lead to a reduction of the effective
        distance. The default value is None, which means it is not
        applied and keep all valid strategies. n_zombie must be set
        to an integer bigger than 0 if used.
        `n_max`: The number of smallest max distance loss values
        being kept. Keeping more than n_max = 1 values helps with
        the global maximization of the effective distance. The default
        value is None, which means it is not applied and we keep all
        valid strategies. n_max must be set to an integer bigger than
        0 if used.
        `n_sum`: The number of smallest sum distance loss values
        being kept. Keeping more than n_sum = 1 values helps with
        the global maximization of the effective distance. The default
         value is None, which means it is not applied and we keep all
        valid strategies. n_sum must be set to an integer bigger than
        0 if used.
        `n_skip`: The minimum number of strategies to evaluate for each
        assignment of link defects to data defects or ancilla zombies
        before moving to the next assignment. We only move on if we
        tried n_skip samples that were sub-optimal and were not recorded.
        The default value is None, which means it is not applied and we
        keep all valid strategies. n_skip must be set to an integer bigger
        than 0 if used.
    """

    link_defect_to_ancilla: bool = False
    """A strong heuristic that restricts the search to the
    strategies that only map link defects to ancilla defects."""
    n_zombie: Optional[int] = None
    """Number of configurations to keep that have the same
    number of zombie qubits, before stopping the search for
    each assignment of link defects.
    """
    n_max: Optional[int] = None
    """Number of minimal max distance loss values to keep. Default
    value is None which corresponds to keeping all."""
    n_sum: Optional[int] = None
    """Number of minimal sum distance loss values to keep. Default
    value is None which corresponds to keeping all."""
    n_skip: Optional[int] = None
    """
    Number of sub-optimal strategies to evaluate for each assignment
    of link defects before moving on to the next assignement. Here
    sub-optimal means that the strategy was not worth recording because
    the distance drop was greater than the ones previously recorded.
    """
```

### TriangularCheck
A class that stores information about an efficient check. The simplest case would have two opposite
triangular checks around a single defective ancilla. More complex scenarios mix these triangular
checks with gauge checks found in the Auger method.

### EfficientDefectsCluster
A class that contains information about a DefectsCluster object on which we applied efficient checks.