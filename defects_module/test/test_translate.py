"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""


def test_memory_circuit(memory_circuit, snapshot):
    """Test the stim output for a logical memory circuit on a nondefective patch."""
    assert memory_circuit.stim_circuit == snapshot
    assert memory_circuit.observables == snapshot
