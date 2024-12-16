"""
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.  

SPDX-License-Identifier: CC-BY-NC-4.0
"""

import warnings
from functools import singledispatch
from itertools import product
from random import shuffle
from typing import Iterator

import numpy
import plotly.graph_objects as go
import scipy

from defects_module.base import PauliT, Pos, Stabilizer, SuperStabilizer
from defects_module.code_library import RotatedSurfaceCode
from defects_module.defects import DefectiveSurfaceCode

# define a color palette for x-type gauges
# in the superstabilizers
x_base_color = "#CE2D40", "#E27985"
x_color_pairs = [
    ("#BD5E49", "#D49587"),
    ("#E55724", "#EE8F6D"),
    ("#DC8118", "#F1BB7E"),
    ("#A5803B", "#C9A769"),
]

# define a color paletter for z-type gauges
# in the superstabilizers
z_base_color = "#128AA5", "#78DAF0"
z_color_pairs = [
    ("#207DD5", "#7BB5EB"),
    ("#564EF9", "#9C97FB"),
    ("#A352E0", "#D2AAF0"),
    ("#5D7598", "#B3BFD0"),
]


def next_color(color_pairs) -> Iterator[tuple[str, str]]:
    """Line and fill color choice pairing.

    We cycle through this list because when care to differentiate amongst nearby components.
    """
    shuffle(color_pairs)
    while True:
        for colors in color_pairs:
            yield colors


def get_stabilizer_color(stabilizer: Stabilizer) -> tuple[str, str]:
    """Default color choice for a specific stabilizer type e.g. X or Z"""
    if stabilizer.type == PauliT.X:
        return x_base_color
    if stabilizer.type == PauliT.Z:
        return z_base_color

    # Mixed stabilizer
    return "DarkOrange", "OliveDrab"


def plot_stabilizer(
    stab: Stabilizer, fig: go.Figure, line_color: str, fillcolor: str, showlegend=False
):
    """Plot the stabilizer with the given line color and fill."""
    ancilla = stab.only_ancilla
    if not isinstance(ancilla, Pos) or not all(
        isinstance(qubit, Pos) for qubit in stab.data_qubits
    ):
        warnings.warn(
            f"Skipping plotting stabilizer {stab} because the qubit's do not "
            "have an associated position.",
            stacklevel=2,
        )
        return

    # For now the toggling only works for unmixed superstabilizers
    if stab.type is PauliT.X:
        legendgroup = "S_X"
    elif stab.type is PauliT.Z:
        legendgroup = "S_Z"
    else:
        legendgroup = None

    # Plot Stabilizer patch
    if stab.weight == 1:
        # Weight one stabilizer is a line connecting from ancilla to data qubits
        xs = [ancilla.x, stab.data_qubits[0].x]
        ys = [ancilla.y, stab.data_qubits[0].y]
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                line=dict(
                    color=line_color,
                    width=2,
                ),
                showlegend=showlegend,
                fillcolor=fillcolor,
                name=legendgroup,
                hoverinfo="none",
                legendgroup=legendgroup,
            )
        )
    elif stab.weight == 2:
        # Weight two stabilizer is a polygon connecting from ancilla to data qubits
        xs = [ancilla.x]
        ys = [ancilla.y]
        xs += [qubit.x for qubit in stab.data_qubits]
        ys += [qubit.y for qubit in stab.data_qubits]
        xs.append(ancilla.x)
        ys.append(ancilla.y)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                line=dict(
                    color=line_color,
                    width=2,
                ),
                showlegend=showlegend,
                fillcolor=fillcolor,
                name=legendgroup,
                hoverinfo="none",
                legendgroup=legendgroup,
            )
        )
    elif stab.weight == 3:
        # Weight three stabilizer is a triangle
        xs = [qubit.x for qubit in stab.data_qubits]
        ys = [qubit.y for qubit in stab.data_qubits]
        xs.append(xs[0])
        ys.append(ys[0])

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                line=dict(
                    color=line_color,
                    width=2,
                ),
                showlegend=showlegend,
                fillcolor=fillcolor,
                name=legendgroup,
                hoverinfo="none",
                legendgroup=legendgroup,
            )
        )
    elif stab.weight == 4:
        # Weight four stabilizer is a rectangle connecting from data qubits to ancilla
        xs = [ancilla.x - 1, ancilla.x + 1, ancilla.x + 1, ancilla.x - 1, ancilla.x - 1]
        ys = [ancilla.y - 1, ancilla.y - 1, ancilla.y + 1, ancilla.y + 1, ancilla.y - 1]

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                line=dict(
                    color=line_color,
                    width=2,
                ),
                showlegend=showlegend,
                fillcolor=fillcolor,
                name=legendgroup,
                hoverinfo="none",
                legendgroup=legendgroup,
            )
        )
    else:
        warnings.warn(
            f"Cannot plot stabilizer {stab} of weight {stab.weight}. "
            "Currently only support weight 2, 3, and 4 stabilizers",
            stacklevel=2,
        )

    # Add data qubit markers
    xs = [qubit.x for qubit in stab.data_qubits]
    ys = [qubit.y for qubit in stab.data_qubits]
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            showlegend=False,
            mode="markers",
            marker_size=8,
            marker_color="#474645",
            name="data",
        )
    )


def stabilizer_center(stab: SuperStabilizer) -> numpy.ndarray:
    """Compute central position of the superstabilizer"""
    # get the positions of all ancillas in the superstabilizer
    ancillas = numpy.array([[ancilla.x, ancilla.y] for ancilla in stab.ancilla])
    # return the average position (i.e. center of mass)
    return numpy.mean(ancillas, axis=0)


def make_clusters(
    superstabilizers: list[SuperStabilizer],
) -> dict[int, list[SuperStabilizer]]:
    """Compute clusters of superstabilizers to assign colors better when plotting"""
    # if no superstabilizers then return an empty dict
    if len(superstabilizers) == 0:
        return dict()
    # get the center of mass for each stabilizer
    pos = numpy.array([stabilizer_center(stab) for stab in superstabilizers])
    # compute the distance matrix needed to compute the linkage matrix
    pdist = scipy.spatial.distance.pdist(pos)
    # compute the linkage matrix
    linkage_matrix = scipy.cluster.hierarchy.ward(pdist)
    # cluster the superstabilizers with the linkage matrix based on cophenetic distance
    # a max distance of 16, corresponding to a cluster with 8 data qubits on its diameter
    # distance could be changed, it is empirical
    indices = scipy.cluster.hierarchy.fcluster(linkage_matrix, t=8, criterion="distance")
    clusters: dict[int, list[SuperStabilizer]] = {j: [] for j in indices}
    for i, j in enumerate(indices):
        clusters[j].append(superstabilizers[i])
    return clusters


def plot_patch(patch: RotatedSurfaceCode, fig: go.Figure):
    """Plot a specific RotatedSurfaceCode patch
    XXX: We have not plotted the logicals. Add this in later.
    """
    # data needed to plot the stabilizers
    stabilizer_data: list[tuple[Stabilizer, str, str]] = []
    # dict where keys are ancillas and values are all checks
    # (undamaged stabilizers or gauge checks) in which they appear
    ancillas: dict[Pos, list[str]] = dict()

    # add undamaged stabilizers to plotting data
    for stab in patch.stabilizers:
        if not isinstance(stab, SuperStabilizer):
            line_color, fillcolor = get_stabilizer_color(stab)
            stabilizer_data.append((stab, line_color, fillcolor))
            ancillas[stab.only_ancilla] = [str(stab.pauli_str)]

    # create clusters for the superstabilizers based on distance
    # such that colors don't overlap when plotting
    clusters = make_clusters(
        [stab for stab in patch.stabilizers if isinstance(stab, SuperStabilizer)]
    )

    # add superstabilizers to plotting data
    for _, superstabilizers in clusters.items():
        # colors generators for x- and z-type
        # superstabilizers
        x_color_gen = next_color(x_color_pairs)
        z_color_gen = next_color(z_color_pairs)
        for stab in superstabilizers:
            if stab.type == PauliT.X:
                line_color, fillcolor = next(x_color_gen)
            else:
                line_color, fillcolor = next(z_color_gen)
            for gauge in stab.gauges:
                stabilizer_data.append((gauge, line_color, fillcolor))
                if gauge.only_ancilla in ancillas.keys():
                    ancillas[gauge.only_ancilla].append(str(gauge.pauli_str))
                else:
                    ancillas[gauge.only_ancilla] = [str(gauge.pauli_str)]

    # sort the plotting data per stabilizer wieght so lower weights stabilizers
    # are on top (i.e. no hidden gauge check)
    stabilizer_data = sorted(stabilizer_data, key=lambda x: len(x[0].data_qubits), reverse=True)

    # plot stabilizers, ordered by weight
    legendshowed = set()
    for gauge_or_stab, line_color, fillcolor in stabilizer_data:
        showlegend = gauge_or_stab.type not in legendshowed
        plot_stabilizer(gauge_or_stab, fig, line_color, fillcolor, showlegend=showlegend)
        if showlegend:
            legendshowed.add(gauge_or_stab.type)

    # Add ancilla qubits markers
    for ancilla, pauli_strs in ancillas.items():
        pauli_str = ",".join(pauli_strs)
        fig.add_trace(
            go.Scatter(
                x=[ancilla.x],
                y=[ancilla.y],
                showlegend=False,
                mode="markers",
                marker_size=10,
                marker_color="#787776",
                name=f"Ancilla {pauli_str}",
                hoverlabel=dict(namelength=-1),
            )
        )

    if isinstance(patch, DefectiveSurfaceCode):
        # Overlay the deformed boundaries
        pauli_types = (
            [PauliT.X, PauliT.Z] if patch.vertical_logical == PauliT.X else [PauliT.Z, PauliT.X]
        )
        line_colors = (
            [x_base_color[0], z_base_color[0]]
            if patch.vertical_logical == PauliT.X
            else [z_base_color[0], x_base_color[0]]
        )
        for label, logical in patch.deformed_logicals.items():
            xs = []
            ys = []
            for q in logical:
                xs.append(q.x)
                ys.append(q.y)
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    showlegend=label in ["right", "top"],
                    line_color=line_colors[label in ["bottom", "top"]],
                    line_width=4,
                    name=f"{pauli_types[label in ['bottom', 'top']]}_L",
                    legendgroup=f"{pauli_types[label in ['bottom', 'top']]}_L",
                )
            )
    else:
        # Overlay the logicals
        for pauli_type, logical in patch.logicals.items():
            xs = []
            ys = []
            for q in logical:
                xs.append(q.x)
                ys.append(q.y)
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    line_width=8,
                    name=f"{pauli_type}_L",
                )
            )


def plot_defects(patch: DefectiveSurfaceCode, fig: go.Figure):
    """Plot the dead qubits and links for an AugerRotatedSurfaceCode."""
    xs = []
    ys = []
    for qubit in patch.qubit_defects:
        xs.append(qubit.x)
        ys.append(qubit.y)

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers",
            showlegend=False,
            name="qubit defect",
            marker=dict(symbol="x", color="red", size=12),
        )
    )

    for link in patch.link_defects:
        fig.add_trace(
            go.Scatter(
                x=[link[0].x, (link[0].x + link[1].x) / 2, link[1].x],
                y=[link[0].y, (link[0].y + link[1].y) / 2, link[1].y],
                mode="lines",
                showlegend=False,
                name="link defect",
                line=dict(color="red", dash="dot", width=2),
            )
        )


def plot_extent(*patches: RotatedSurfaceCode, fig: go.Figure):
    """Plot the spatial extent"""
    xmin, xmax, ymin, ymax = get_extent(*patches)
    xs_data: list[int] = []
    xs_ancilla: list[int] = []
    ys_data: list[int] = []
    ys_ancilla: list[int] = []

    for x, y in product(range(xmin, xmax + 1, 2), range(ymin, ymax + 1, 2)):
        xs_ancilla.append(x)
        ys_ancilla.append(y)

    for x, y in product(range(xmin + 1, xmax, 2), range(ymin + 1, ymax, 2)):
        xs_data.append(x)
        ys_data.append(y)

    fig.add_trace(
        go.Scatter(
            x=xs_data,
            y=ys_data,
            showlegend=False,
            mode="markers",
            marker_size=8,
            marker_color="#474645",
            name="qubit",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs_ancilla,
            y=ys_ancilla,
            showlegend=False,
            mode="markers",
            marker_size=10,
            marker_color="#787776",
            name="qubit",
        )
    )


def get_extent(*patches: RotatedSurfaceCode) -> tuple[int, int, int, int]:
    """Get the total spatial extent of a set of patches"""
    min_x = 0
    max_x = 0
    min_y = 0
    max_y = 0
    for patch in patches:
        min_x = min(min_x, patch.extent[0])
        max_x = max(max_x, patch.extent[1])
        min_y = min(min_y, patch.extent[2])
        max_y = max(max_y, patch.extent[3])

    return min_x, max_x, min_y, max_y


@singledispatch
def plot(obj):
    """Default implementaion."""
    raise NotImplementedError(f"Cannot plot {obj}")


@plot.register
def _(*patches: RotatedSurfaceCode, style: str = "both"):
    """Plot a collection of RotatedSurfaceCode patches, including data qubits, ancillas, and
    logicals.
    """
    assert style in [
        "extent",
        "patch",
        "both",
    ], "Only 'extent', 'patch', and 'both' arguments supported for style"
    if not all(isinstance(patch, RotatedSurfaceCode) for patch in patches):
        raise NotImplementedError("Only RotatedSquare surface codes currently supported.")

    fig = go.Figure()

    if style in ["extent", "both"]:
        plot_extent(*patches, fig=fig)

    if style in ["patch", "both"]:
        for patch in patches:
            plot_patch(patch, fig)
            if isinstance(patch, DefectiveSurfaceCode):
                plot_defects(patch, fig)

    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False, scaleanchor="x")
    fig.update_layout(autosize=False, width=800)

    return fig
