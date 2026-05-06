"""
Relevant paper on Chaos Game Representation in Bioinformatics
https://www.sciencedirect.com/science/article/pii/S2001037021004736
"""

import math
import numpy as np
import numpy.typing as npt


def generate_square_points(
    num_vertices: int,
    radius: float = 1.0,
    center: tuple[float, float] = (0.0, 0.0),
    rotation_deg: int = 0,
) -> list[tuple[float, float]]:
    cx, cy = center
    points = []

    for i in range(num_vertices):
        x = radius * math.sin(
            (2 * math.pi * i / num_vertices) + math.radians(rotation_deg)
        )
        y = radius * math.cos(
            (2 * math.pi * i / num_vertices) + math.radians(rotation_deg)
        )
        points.append(
            (round(x, 6), round(y, 6))
        )  # Rounding to prevent super small fractions

    return points


def label_chaos_points(
    vertex_names: str | list[str], vertices: list[tuple[float, float]]
) -> dict[str, tuple[float, float]] | None:
    if len(vertex_names) != len(vertices):
        return
    uppercase = [label.upper() for label in vertex_names]
    lowercase = [label.lower() for label in vertex_names]
    case_insensitive = uppercase + lowercase
    doubled_vertices = vertices * 2
    return {
        label: vertex for (label, vertex) in zip(case_insensitive, doubled_vertices)
    }


def create_chaos(
    sequence: str,
    center: tuple[float, float],
    chaos_dict: dict[str, tuple[float, float]],
) -> list[tuple[float, float]] | None:
    num_vertices = len(chaos_dict)
    m = num_vertices // 4
    if num_vertices == 4:
        scaling_factor = 0.5
    else:
        scaling_factor = 1 - (
            math.sin(math.pi / num_vertices)
            / (
                math.sin(math.pi / num_vertices)
                + math.sin(math.pi / num_vertices + (2 * math.pi * m / num_vertices))
            )
        )
    cgr = []
    cgr_marker = center[:]
    for letter in sequence:
        step_direction = chaos_dict[letter]
        if step_direction:
            cgr_marker = (
                (
                    cgr_marker[0]
                    + (scaling_factor * (step_direction[0] - cgr_marker[0]))
                ),
                (
                    cgr_marker[1]
                    + (scaling_factor * (step_direction[1] - cgr_marker[1]))
                ),
            )
            cgr.append(cgr_marker)
        else:
            return None
    return cgr


def cgr_to_fcgr(
    cgr: list[tuple[float, float]], resolution: int, radius: float
) -> npt.NDArray[np.intc]:
    """
    Converts coordinates of Chaos Game Representation to numpy grid as
    Frequency Chaos Game Representation. Allows for fixed size for NN input.
    """
    grid = np.zeros((resolution, resolution))
    coords = np.array(cgr)
    x_min, x_max = -radius, radius
    y_min, y_max = -radius, radius

    for x, y in coords:
        col = int((x - x_min) / (x_max - x_min) * (resolution - 1))
        row = int((y - y_min) / (y_max - y_min) * (resolution - 1))
        col = max(0, min(resolution - 1, col))
        row = max(0, min(resolution - 1, row))
        grid[row, col] += 1
    return grid
