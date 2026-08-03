"""Six-branch electrostatic dispersion solvers.

This module contains the reduced cold three-species polynomial solver and the
cold-ion/two-electron full-Maxwellian solver. Both solvers return
``SixComplexRoots`` and preserve complete-complex branch continuity on the two
signed wavenumber half-axes. The module also provides the HDF5 exporter used by
the plotting driver.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, permutations
from pathlib import Path
from typing import Callable, Iterable
import warnings

import h5py
import numpy
from scipy.optimize import root
from scipy.special import wofz


############################################################
# Constants and branch conventions
############################################################

CGS = {
    "me": 9.1094e-28,
    "mp": 1.6726231e-24,
    "mpi0": 2.406175207e-25,
    "e": 4.8032068e-10,
    "c": 29979245800.0,
}

BRANCH_LABELS = (
    r"$\omega_1$: positive-frequency outer sheet",
    r"$\omega_2$: negative-frequency outer sheet",
    r"$\omega_3$: $G_2$ first continuous sheet",
    r"$\omega_4$: $G_2$ second continuous sheet",
    r"$\omega_5$: $G_3$ first continuous sheet",
    r"$\omega_6$: $G_3$ second continuous sheet",
)

GROUP_NAMES = ("G1_outer", "G2_positive_drift", "G3_negative_drift")
GROUP_COLUMNS = ((0, 1), (2, 3), (4, 5))

# Zero-based source branch at +k for each branch at -k:
# roots(-k, j) = -conj(roots(+k, OPPOSITE_K_PARTNER[j])).
OPPOSITE_K_PARTNER = numpy.asarray([1, 0, 2, 3, 4, 5], dtype=int)


__all__ = [
    "BRANCH_LABELS",
    "CGS",
    "OPPOSITE_K_PARTNER",
    "SixComplexRoots",
    "save_solution_hdf5",
    "solve_six_full_maxwellian_branches",
    "solve_six_reduced_branches",
]

# Dimensionless tolerances after normalization by the total plasma frequency.
_NONREAL_TOL = 2.0e-10
_CANONICAL_REAL_ORDER_TOL = 2.0e-10
_PAIR_CONJUGACY_TOL = 2.0e-7
_PAIR_ASSIGNMENT_TIE_ATOL = 2.0e-7
_PAIR_ASSIGNMENT_TIE_RTOL = 2.0e-6
_GROWTH_COMPARISON_TOL = 2.0e-10

# Three partitions of four inner roots into two unordered pairs.
_PAIR_PARTITIONS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


@dataclass(frozen=True)
class HalfAxisSolution:
    """Direct roots on one open half-axis, ordered by increasing |k|."""

    k_path: numpy.ndarray
    roots: numpy.ndarray
    raw_roots: numpy.ndarray
    raw_source_index: numpy.ndarray
    conjugate_partition_enforced: numpy.ndarray
    analytic_tie_break_used: numpy.ndarray


@dataclass(frozen=True)
class SixComplexRoots:
    """Six complete-complex roots and solver diagnostics.

    The same result class is used by the reduced-polynomial and full-Maxwellian
    solvers.  ``k_bs`` has shape (2, 3): rows correspond to the two inner pairs
    and columns correspond to the negative-side, positive-side, and common
    branch boundaries.
    """

    k: numpy.ndarray
    roots: numpy.ndarray
    k_bs: numpy.ndarray
    parameters: dict[str, object] = field(default_factory=dict)

    raw_roots: numpy.ndarray | None = None
    raw_source_index: numpy.ndarray | None = None
    reduced_seed: numpy.ndarray | None = None
    solved_directly: numpy.ndarray | None = None
    relative_residual: numpy.ndarray | None = None
    converged: numpy.ndarray | None = None
    iterations: numpy.ndarray | None = None
    candidate_count: numpy.ndarray | None = None
    pseudoarclength_used: numpy.ndarray | None = None
    contour_used: numpy.ndarray | None = None
    contour_certified: numpy.ndarray | None = None
    contour_count: numpy.ndarray | None = None
    analytic_tie_break_used: numpy.ndarray | None = None
    conjugate_partition_enforced: numpy.ndarray | None = None

    positive_anchor_index: int | None = None
    negative_anchor_index: int | None = None
    positive_homotopy: object | None = None
    negative_homotopy: object | None = None

    same_real_error: numpy.ndarray | None = None
    opposite_imag_error: numpy.ndarray | None = None
    growth_magnitude_error: numpy.ndarray | None = None
    complex_pair_mask: numpy.ndarray | None = None
    pair_condition_mask: numpy.ndarray | None = None

    opposite_k_conversion_error: numpy.ndarray | None = None
    opposite_k_root_set_mismatch: float = numpy.nan
    complete_root_assignment_error: float = numpy.nan

    k_b1_negative: float = numpy.nan
    k_b1_positive: float = numpy.nan
    k_b1_common: float = numpy.nan
    k_b2_negative: float = numpy.nan
    k_b2_positive: float = numpy.nan
    k_b2_common: float = numpy.nan
    reference_k_bs: numpy.ndarray | None = None

    full_k_b1_negative: float = numpy.nan
    full_k_b1_positive: float = numpy.nan
    full_k_b1_common: float = numpy.nan
    full_k_b2_negative: float = numpy.nan
    full_k_b2_positive: float = numpy.nan
    full_k_b2_common: float = numpy.nan

    growth_hierarchy_mask: numpy.ndarray | None = None
    growth_hierarchy_margin: numpy.ndarray | None = None
    minimum_growth_hierarchy_margin: float = numpy.nan
    growth_hierarchy_violation_count: int = 0

    slope_boundary_k: numpy.ndarray | None = None
    slope_orders: numpy.ndarray | None = None
    slope_real_left: numpy.ndarray | None = None
    slope_real_right: numpy.ndarray | None = None
    slope_imag_left: numpy.ndarray | None = None
    slope_imag_right: numpy.ndarray | None = None
    slope_real_relative_mismatch: numpy.ndarray | None = None
    slope_imag_relative_mismatch: numpy.ndarray | None = None
    maximum_rule5_relative_mismatch: float = numpy.nan

    def __post_init__(self) -> None:
        if self.k.ndim != 1:
            raise ValueError("k must be one-dimensional.")
        if self.roots.shape != (self.k.size, 6):
            raise ValueError("roots must have shape (len(k), 6).")
        if self.k_bs.shape != (2, 3):
            raise ValueError("k_bs must have shape (2, 3).")

    @property
    def roots_real(self) -> numpy.ndarray:
        return self.roots.real

    @property
    def roots_imag(self) -> numpy.ndarray:
        return self.roots.imag

    @property
    def k_b_common(self) -> float:
        finite = self.k_bs[:, 2]
        finite = finite[numpy.isfinite(finite)]
        if finite.size == 0:
            return numpy.nan
        return float(numpy.min(finite))


############################################################
# Polynomial construction and direct root solution
############################################################


def wave_equation_two_electrons_ions_coefficients(
    k: float | complex,
    wpe_negative: float,
    vd_negative: float,
    wpe_positive: float,
    vd_positive: float,
    wpi: float,
) -> numpy.ndarray:
    """Return sixth-order coefficients in descending powers of omega."""

    a6 = 1.0

    a5 = -2.0 * k * (vd_negative + vd_positive)

    a4 = (
        k * k * ( vd_negative * vd_negative + 4.0 * vd_negative * vd_positive + vd_positive * vd_positive)
        - wpi * wpi
        - wpe_negative * wpe_negative
        - wpe_positive * wpe_positive
    )

    a3 = 2.0 * k * (
        wpe_negative * wpe_negative * vd_positive
        + wpe_positive * wpe_positive * vd_negative
        + wpi * wpi * (vd_negative + vd_positive)
        - k * k * vd_negative * vd_positive * (vd_negative + vd_positive)
    )

    a2 = (
        k**4 * vd_negative**2 * vd_positive**2
        - k * k * wpe_negative**2 * vd_positive**2
        - k * k * wpe_positive**2 * vd_negative**2
        - k * k * wpi * wpi * (vd_negative**2 + 4.0 * vd_negative * vd_positive + vd_positive**2)
    )

    a1 = 2.0 * k**3 * wpi * wpi * vd_negative * vd_positive * (vd_negative + vd_positive)

    a0 = -k**4 * wpi**2 * vd_negative**2 * vd_positive**2

    dtype = numpy.complex128 if numpy.iscomplexobj(k) else float

    return numpy.asarray([a6, a5, a4, a3, a2, a1, a0], dtype=dtype)


def _frequency_scale(
    wpe_negative: float,
    wpe_positive: float,
    wpi: float,
) -> float:
    scale = float(numpy.sqrt(wpe_negative * wpe_negative + wpe_positive * wpe_positive + wpi * wpi))
    if not numpy.isfinite(scale) or scale <= 0.0:
        raise ValueError("At least one plasma frequency must be positive.")
    return scale


def _validate_species_order(vd_negative: float, vd_positive: float) -> None:
    if not vd_negative < 0.0:
        raise ValueError("vd_negative must be strictly negative.")
    if not vd_positive > 0.0:
        raise ValueError("vd_positive must be strictly positive.")


def reduced_roots_at_k(
    k: float,
    wpe_negative: float,
    vd_negative: float,
    wpe_positive: float,
    vd_positive: float,
    wpi: float,
) -> numpy.ndarray:
    """Directly solve all six unordered roots at one nonzero real k.

    The polynomial is solved in x=omega/omega_scale. This scaling improves
    coefficient conditioning without changing any physical root.
    """

    if k == 0.0:
        raise ValueError("Use reduced_limiting_roots_at_zero at k=0.")

    scale = _frequency_scale(wpe_negative, wpe_positive, wpi)
    coefficients = wave_equation_two_electrons_ions_coefficients(k, wpe_negative, vd_negative, wpe_positive, vd_positive, wpi,).astype(numpy.complex128)

    # If omega=scale*x, division by scale**6 gives coefficient i divided by
    # scale**i for coefficients ordered from omega**6 to omega**0.
    scaled_coefficients = coefficients / scale ** numpy.arange(7, dtype=float)
    roots = numpy.roots(scaled_coefficients).astype(numpy.complex128) * scale

    if roots.shape != (6,):
        raise RuntimeError(f"Expected six roots at k={k!r}; got {roots.shape}.")
    if not (numpy.isfinite(roots.real).all() and numpy.isfinite(roots.imag).all()):
        raise RuntimeError(f"Nonfinite polynomial root at k={k!r}.")
    return roots


def reduced_limiting_roots_at_zero(
    wpe_negative: float,
    wpe_positive: float,
    wpi: float,
) -> numpy.ndarray:
    """Return exact k->0 limiting roots in the prescribed branch order."""

    total = _frequency_scale(wpe_negative, wpe_positive, wpi)

    return numpy.asarray([+total, -total, 0.0j, 0.0j, 0.0j, 0.0j], dtype=numpy.complex128,)


############################################################
# Complete-complex pair continuation
############################################################


def _reduced_pair_descriptor(pair: numpy.ndarray) -> tuple[complex, complex, float]:
    """Return pair center, squared half-separation, and RMS radius.

    The squared half-separation is analytic through a conjugate-to-real
    coalescence, whereas the two individual square-root sheets are not.
    """

    center = 0.5 * (pair[0] + pair[1])
    half_difference = 0.5 * (pair[0] - pair[1])
    squared_half_separation = half_difference * half_difference
    radius = float(numpy.sqrt(0.5 * (abs(pair[0] - center) ** 2 + abs(pair[1] - center) ** 2)))

    return center, squared_half_separation, radius


def _secant_predict(
    current: numpy.ndarray,
    previous: numpy.ndarray | None,
    step_ratio: float,
) -> numpy.ndarray:
    if previous is None:
        return current.copy()
    return current + step_ratio * (current - previous)


def _multiscale_secant_predict(
    values: numpy.ndarray,
    grid: numpy.ndarray,
    next_index: int,
    maximum_order: int,
) -> numpy.ndarray:
    """Predict one step using secants with n=1,...,maximum_order.

    The predictor is the inverse-order weighted mean of complete-complex
    extrapolations. It directly encodes the Rule-5 preference that left and
    right discrete slopes agree over several stencil widths. No real/imaginary
    component is predicted or assigned separately.
    """

    current_index = next_index - 1
    current = values[current_index]
    if current_index <= 0 or maximum_order <= 0:
        return current.copy()

    predictions = []
    weights = []
    for order in range(1, min(maximum_order, current_index) + 1):
        earlier_index = current_index - order
        denominator = grid[current_index] - grid[earlier_index]
        if denominator == 0.0:
            continue
        slope = (current - values[earlier_index]) / denominator
        prediction = current + slope * (grid[next_index] - grid[current_index])
        predictions.append(prediction)
        weights.append(1.0 / order)

    if not predictions:
        return current.copy()
    weight_array = numpy.asarray(weights, dtype=float)
    stacked = numpy.stack(predictions, axis=0)
    return numpy.tensordot(weight_array / weight_array.sum(), stacked, axes=(0, 0))


def _select_outer_pair(
    raw_roots: numpy.ndarray,
    predictors: numpy.ndarray,
    scale: float,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Select ordered G1 roots by complete-complex continuity."""

    best: tuple[float, int, int] | None = None
    for first, second in permutations(range(6), 2):
        cost = float((abs(predictors[0] - raw_roots[first]) + abs(predictors[1] - raw_roots[second]))/ scale)

        candidate = (cost, first, second)

        if best is None or candidate < best:
            best = candidate

    if best is None:
        raise RuntimeError("No outer-pair assignment was found.")
    source = numpy.asarray([best[1], best[2]], dtype=int)
    return raw_roots[source], source


def _anchor_outer_pair(
    raw_roots: numpy.ndarray,
    scale: float,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    targets = numpy.asarray([+scale, -scale], dtype=numpy.complex128)
    return _select_outer_pair(raw_roots, targets, scale)


def _candidate_inner_partitions(
    inner_roots: numpy.ndarray,
    scale: float,
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    """Return physically admissible pair partitions of four inner roots.

    Nonreal roots are kept with their complex conjugates. This prevents a
    growing root from being paired with a root on an unrelated analytic sheet.
    """

    nonreal = numpy.abs(inner_roots.imag) > _NONREAL_TOL * scale
    number_nonreal = int(numpy.count_nonzero(nonreal))
    if number_nonreal not in (2, 4):
        return _PAIR_PARTITIONS

    scored: list[tuple[float, tuple[tuple[int, int], tuple[int, int]]]] = []

    for partition in _PAIR_PARTITIONS:
        valid = True
        conjugacy_error = 0.0
        for first, second in partition:
            if nonreal[first] != nonreal[second]:
                valid = False
                break
            if nonreal[first]:
                conjugacy_error += float(abs(inner_roots[first] - numpy.conjugate(inner_roots[second]))/ scale)
        
        if valid:
            scored.append((conjugacy_error, partition))

    if not scored:
        return _PAIR_PARTITIONS

    scored.sort(key=lambda item: (item[0], item[1]))
    best_error = scored[0][0]
    tolerance = _PAIR_CONJUGACY_TOL + 10.0 * best_error
    admissible = tuple(partition for error, partition in scored if error <= tolerance)
    return admissible if admissible else (scored[0][1],)


def _pair_set_cost(
    predicted_pair: numpy.ndarray,
    candidate_pair: numpy.ndarray,
    scale: float,
) -> float:
    """Minimum intact-complex matching cost between two unordered pairs."""

    direct = abs(predicted_pair[0] - candidate_pair[0]) + abs(predicted_pair[1] - candidate_pair[1])

    crossed = abs(predicted_pair[0] - candidate_pair[1]) + abs(predicted_pair[1] - candidate_pair[0])

    return float(min(direct, crossed) / scale)


def _order_self_connected_pair(
    pair: numpy.ndarray,
    source: numpy.ndarray,
    sign: int,
    scale: float,
    predicted_pair: numpy.ndarray | None = None,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Order a G2/G3 pair by complete-root continuity.

    At the one-sided anchor, a nonreal pair is oriented with the growing member
    first.  At subsequent points, the direct and crossed complete-root
    assignments are compared with the multiscale predictor.  Componentwise or
    pointwise |Re(omega)| sorting is not applied after a nondegenerate
    continuity decision.
    """

    nonreal = numpy.max(numpy.abs(pair.imag)) > _NONREAL_TOL * scale
    if predicted_pair is None:
        order = numpy.argsort(-pair.imag) if nonreal else numpy.argsort(sign * pair.real)
        return pair[order], source[order]

    direct = abs(pair[0] - predicted_pair[0]) + abs(pair[1] - predicted_pair[1])
    crossed = abs(pair[1] - predicted_pair[0]) + abs(pair[0] - predicted_pair[1])
    tolerance = (
        _PAIR_ASSIGNMENT_TIE_ATOL
        + _PAIR_ASSIGNMENT_TIE_RTOL * max(1.0, float(min(direct, crossed) / scale))
    ) * scale
    if direct + tolerance < crossed:
        order = numpy.asarray([0, 1], dtype=int)
    elif crossed + tolerance < direct:
        order = numpy.asarray([1, 0], dtype=int)
    elif nonreal:
        order = numpy.argsort(-pair.imag)
    else:
        order = numpy.argsort(sign * pair.real)
    return pair[order], source[order]


def _inner_assignment_options(
    inner_roots: numpy.ndarray,
    inner_source: numpy.ndarray,
    partitions: Iterable[tuple[tuple[int, int], tuple[int, int]]],
    predicted_roots: numpy.ndarray | None,
    predicted_centers: numpy.ndarray | None,
    predicted_squared_separations: numpy.ndarray | None,
    scale: float,
) -> list[dict]:
    options: list[dict] = []
    anchor = predicted_roots is None

    for partition in partitions:
        for swap in (False, True):
            assigned_partition = ((partition[1], partition[0]) if swap else partition)

            descriptors = []
            total_cost = 0.0

            for group_index, local_pair in enumerate(assigned_partition):
                values = inner_roots[list(local_pair)]
                center, squared_separation, radius = _reduced_pair_descriptor(values)
                descriptors.append((center, squared_separation, radius))

                if not anchor:
                    branch_slice = slice(2 * group_index, 2 * group_index + 2)
                    total_cost += _pair_set_cost(predicted_roots[branch_slice], values, scale)
                    total_cost += 0.75 * float(abs(center - predicted_centers[group_index]) / scale)
                    total_cost += 0.25 * float(abs(squared_separation - predicted_squared_separations[group_index])/ (scale * scale))

            options.append(
                {
                    "cost": float(total_cost),
                    "partition": assigned_partition,
                    "descriptors": descriptors,
                    "source": inner_source,
                }
            )
    return options


def _select_inner_groups(
    inner_roots: numpy.ndarray,
    inner_source: numpy.ndarray,
    sign: int,
    scale: float,
    k_value: float,
    vd_negative: float,
    vd_positive: float,
    *,
    predicted_roots: numpy.ndarray | None = None,
    predicted_centers: numpy.ndarray | None = None,
    predicted_squared_separations: numpy.ndarray | None = None,
) -> tuple[numpy.ndarray, numpy.ndarray, bool, bool]:
    """Select G2 then G3 while preserving complete-complex pair continuity.

    At the one-sided near-zero anchor, G2 is the stronger conjugate pair and G3
    is the weaker pair. Drift-line proximity is only a lower-priority tie-break.
    At later points, pair continuity has priority. A deterministic tie-break is
    used only when coalescence/reconnection makes the continuation costs equal.
    """

    anchor = predicted_roots is None
    partitions = _candidate_inner_partitions(inner_roots, scale)
    conjugate_partition_enforced = len(partitions) < len(_PAIR_PARTITIONS)
    options = _inner_assignment_options(
        inner_roots,
        inner_source,
        partitions,
        predicted_roots,
        predicted_centers,
        predicted_squared_separations,
        scale,
    )
    if not options:
        raise RuntimeError("No admissible G2/G3 assignment was found.")

    analytic_tie_break_used = False
    if anchor:
        positive_target = k_value * vd_positive
        negative_target = k_value * vd_negative

        # Lexicographic selection: larger G2 radius first, then place its center
        # nearer the positive-drift line and the G3 center nearer the negative-
        # drift line. The radius criterion resolves the symmetric case.
        selected = max(
            options,
            key=lambda option: (
                (option["descriptors"][0][2] - option["descriptors"][1][2])/ scale,
                -(abs(option["descriptors"][0][0].real - positive_target) + abs(option["descriptors"][1][0].real - negative_target))/ scale,
                sign * (option["descriptors"][0][0].real - option["descriptors"][1][0].real)/ scale,
            ),
        )
    else:
        options.sort(key=lambda option: option["cost"])
        best_cost = options[0]["cost"]
        tolerance = _PAIR_ASSIGNMENT_TIE_ATOL + _PAIR_ASSIGNMENT_TIE_RTOL * max(1.0, best_cost)

        tied = [option for option in options if option["cost"] <= best_cost + tolerance]

        if len(tied) == 1:
            selected = tied[0]
        else:
            # This is used only at a coalescence or analytic reconnection. It
            # does not override a nondegenerate continuity decision.
            analytic_tie_break_used = True
            selected = max(
                tied,
                key=lambda option: (
                    sign * (option["descriptors"][0][0].real - option["descriptors"][1][0].real)/ scale,
                    (option["descriptors"][0][2] - option["descriptors"][1][2])/ scale,
                ),
            )

    assigned = numpy.empty(4, dtype=numpy.complex128)
    assigned_source = numpy.empty(4, dtype=int)
    for group_index, local_pair in enumerate(selected["partition"]):
        values = inner_roots[list(local_pair)]
        source = inner_source[list(local_pair)]
        branch_slice = slice(2 * group_index, 2 * group_index + 2)
        predicted_pair = (None if predicted_roots is None else predicted_roots[branch_slice])
        values, source = _order_self_connected_pair(values, source, sign, scale, predicted_pair)
        assigned[branch_slice] = values
        assigned_source[branch_slice] = source

    return (
        assigned,
        assigned_source,
        conjugate_partition_enforced,
        analytic_tie_break_used,
    )


def solve_reduced_half_axis(
    k_path: numpy.ndarray,
    wpe_negative: float,
    vd_negative: float,
    wpe_positive: float,
    vd_positive: float,
    wpi: float,
    *,
    sign: int,
    slope_order_count: int = 4,
) -> HalfAxisSolution:
    """Solve one open half-axis from small to large |k|.

    For sign=+1, k_path must be +dk,+2dk,... .
    For sign=-1, k_path must be -dk,-2dk,... .

    The two half-axes are independent. Opposite-k roots are not read or used.
    """

    k_path = numpy.asarray(k_path, dtype=float)
    if k_path.ndim != 1 or k_path.size == 0:
        raise ValueError("k_path must be a nonempty one-dimensional array.")
    if sign not in (-1, +1):
        raise ValueError("sign must be +1 or -1.")
    if slope_order_count < 1:
        raise ValueError("slope_order_count must be at least one.")
    if not numpy.all(sign * k_path > 0.0):
        raise ValueError("k_path contains the wrong sign.")
    if numpy.any(numpy.diff(numpy.abs(k_path)) <= 0.0):
        raise ValueError("k_path must be ordered by increasing |k|.")

    n_k = k_path.size
    scale = _frequency_scale(wpe_negative, wpe_positive, wpi)
    roots = numpy.empty((n_k, 6), dtype=numpy.complex128)
    raw_roots = numpy.empty_like(roots)
    raw_source_index = numpy.empty((n_k, 6), dtype=int)
    conjugate_enforced = numpy.zeros(n_k, dtype=bool)
    tie_used = numpy.zeros(n_k, dtype=bool)

    pair_centers = numpy.empty((n_k, 2), dtype=numpy.complex128)
    pair_squared_separations = numpy.empty((n_k, 2), dtype=numpy.complex128)

    # One-sided near-zero anchor. The root set itself is solved directly at the
    # first nonzero k; only its group names are initialized here.
    raw = reduced_roots_at_k(
        k_path[0],
        wpe_negative,
        vd_negative,
        wpe_positive,
        vd_positive,
        wpi,
    )
    raw_roots[0] = raw
    outer, outer_source = _anchor_outer_pair(raw, scale)
    roots[0, 0:2] = outer
    raw_source_index[0, 0:2] = outer_source

    outer_set = set(int(value) for value in outer_source)
    remaining_source = numpy.asarray([index for index in range(6) if index not in outer_set], dtype=int)
    inner, inner_source, enforced, used_tie = _select_inner_groups(
        raw[remaining_source],
        remaining_source,
        sign,
        scale,
        k_path[0],
        vd_negative,
        vd_positive,
    )
    roots[0, 2:6] = inner
    raw_source_index[0, 2:6] = inner_source
    conjugate_enforced[0] = enforced
    tie_used[0] = used_tie

    for group_index in range(2):
        branch_slice = slice(2 + 2 * group_index, 4 + 2 * group_index)
        (
            pair_centers[0, group_index],
            pair_squared_separations[0, group_index],
            _,
        ) = _reduced_pair_descriptor(roots[0, branch_slice])

    for index in range(1, n_k):
        raw = reduced_roots_at_k(
            k_path[index],
            wpe_negative,
            vd_negative,
            wpe_positive,
            vd_positive,
            wpi,
        )
        raw_roots[index] = raw

        # Rule 5: use several complete-complex secant stencils rather than a
        # special-case smoother at a detected breakpoint.
        predicted_roots = _multiscale_secant_predict(roots, k_path, index, slope_order_count)
        predicted_centers = _multiscale_secant_predict(pair_centers, k_path, index, slope_order_count)
        predicted_squared_separations = _multiscale_secant_predict(pair_squared_separations, k_path, index, slope_order_count)

        outer, outer_source = _select_outer_pair(raw, predicted_roots[0:2], scale)

        roots[index, 0:2] = outer
        raw_source_index[index, 0:2] = outer_source

        outer_set = set(int(value) for value in outer_source)
        remaining_source = numpy.asarray([raw_index for raw_index in range(6) if raw_index not in outer_set], dtype=int,)
        inner, inner_source, enforced, used_tie = _select_inner_groups(
            raw[remaining_source],
            remaining_source,
            sign,
            scale,
            k_path[index],
            vd_negative,
            vd_positive,
            predicted_roots=predicted_roots[2:6],
            predicted_centers=predicted_centers,
            predicted_squared_separations=predicted_squared_separations,
        )
        roots[index, 2:6] = inner
        raw_source_index[index, 2:6] = inner_source
        conjugate_enforced[index] = enforced
        tie_used[index] = used_tie

        for group_index in range(2):
            branch_slice = slice(2 + 2 * group_index, 4 + 2 * group_index)
            (
                pair_centers[index, group_index],
                pair_squared_separations[index, group_index],
                _,
            ) = _reduced_pair_descriptor(roots[index, branch_slice])

    return HalfAxisSolution(
        k_path=k_path.copy(),
        roots=roots,
        raw_roots=raw_roots,
        raw_source_index=raw_source_index,
        conjugate_partition_enforced=conjugate_enforced,
        analytic_tie_break_used=tie_used,
    )


############################################################
# Diagnostics and full signed-grid assembly
############################################################


def reduced_relative_residual(
    k: float,
    omega: complex,
    wpe_negative: float,
    vd_negative: float,
    wpe_positive: float,
    vd_positive: float,
    wpi: float,
) -> float:
    """Return a cancellation-aware residual of the rational equation."""

    if k == 0.0:
        return numpy.nan
    denominators = (
        omega,
        omega - k * vd_negative,
        omega - k * vd_positive,
    )
    if any(abs(value) == 0.0 for value in denominators):
        return numpy.inf

    terms = (
        wpi * wpi / (omega * omega),
        wpe_negative * wpe_negative / (omega - k * vd_negative) ** 2,
        wpe_positive * wpe_positive / (omega - k * vd_positive) ** 2,
    )
    value = 1.0 - terms[0] - terms[1] - terms[2]
    return float(abs(value) / (1.0 + sum(abs(term) for term in terms)))


def _pair_consistency(
    roots: numpy.ndarray,
    scale: float,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    same_real = numpy.empty((roots.shape[0], 2), dtype=float)
    opposite_imag = numpy.empty_like(same_real)
    complex_mask = numpy.empty_like(same_real, dtype=bool)

    for group_index, first_branch in enumerate((2, 4)):
        first = roots[:, first_branch]
        second = roots[:, first_branch + 1]
        same_real[:, group_index] = numpy.abs(first.real - second.real) / scale
        opposite_imag[:, group_index] = numpy.abs(first.imag + second.imag) / scale
        complex_mask[:, group_index] = (numpy.maximum(numpy.abs(first.imag), numpy.abs(second.imag)) > _NONREAL_TOL * scale) & (same_real[:, group_index] <= _PAIR_CONJUGACY_TOL) & (opposite_imag[:, group_index] <= _PAIR_CONJUGACY_TOL)

    return same_real, opposite_imag, complex_mask


def _contiguous_complex_boundary(
    k_side: numpy.ndarray,
    complex_mask: numpy.ndarray,
) -> float:
    """Estimate the outer boundary of a complex regime contiguous with k=0."""

    order = numpy.argsort(numpy.abs(k_side))
    magnitudes = numpy.abs(k_side[order])
    mask = complex_mask[order]
    if magnitudes.size == 0 or not mask[0]:
        return 0.0

    false_indices = numpy.flatnonzero(~mask)
    if false_indices.size == 0:
        # The pair is still complex at the edge of the requested interval, so
        # no physical coalescence boundary has been detected.
        return numpy.nan

    first_false = int(false_indices[0])
    if first_false == 0:
        return 0.0
    return float(0.5 * (magnitudes[first_false - 1] + magnitudes[first_false]))


def _reduced_rule5_slope_diagnostics(
    k: numpy.ndarray,
    roots: numpy.ndarray,
    boundary_k: numpy.ndarray,
    slope_order_count: int,
) -> tuple[
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    float,
]:
    """Evaluate left/right slopes at -kb1,+kb1,-kb2,+kb2.

    For a boundary not reached within the numerical interval, its row remains
    NaN. The relative mismatch uses a velocity floor of 1e-12*c. This function
    is diagnostic only and never changes a polynomial root.
    """

    orders = numpy.arange(1, slope_order_count + 1, dtype=int)
    shape = (4, slope_order_count, 6)
    real_left = numpy.full(shape, numpy.nan, dtype=float)
    real_right = numpy.full(shape, numpy.nan, dtype=float)
    imag_left = numpy.full(shape, numpy.nan, dtype=float)
    imag_right = numpy.full(shape, numpy.nan, dtype=float)
    real_mismatch = numpy.full(shape, numpy.nan, dtype=float)
    imag_mismatch = numpy.full(shape, numpy.nan, dtype=float)
    velocity_floor = 1.0e-12 * CGS["c"]

    for boundary_index, requested_k in enumerate(boundary_k):
        if not numpy.isfinite(requested_k) or requested_k == 0.0:
            continue
        center = int(numpy.argmin(numpy.abs(k - requested_k)))
        for order_index, order in enumerate(orders):
            left_index = center - int(order)
            right_index = center + int(order)
            if left_index < 0 or right_index >= k.size:
                continue
            left_denominator = k[center] - k[left_index]
            right_denominator = k[right_index] - k[center]
            if left_denominator == 0.0 or right_denominator == 0.0:
                continue

            left = (roots[center] - roots[left_index]) / left_denominator
            right = (roots[right_index] - roots[center]) / right_denominator
            real_left[boundary_index, order_index] = left.real
            real_right[boundary_index, order_index] = right.real
            imag_left[boundary_index, order_index] = left.imag
            imag_right[boundary_index, order_index] = right.imag

            real_scale = numpy.maximum.reduce((numpy.abs(left.real), numpy.abs(right.real), numpy.full(6, velocity_floor),))
            imag_scale = numpy.maximum.reduce((numpy.abs(left.imag), numpy.abs(right.imag), numpy.full(6, velocity_floor),))
            real_mismatch[boundary_index, order_index] = (numpy.abs(left.real - right.real) / real_scale)
            imag_mismatch[boundary_index, order_index] = (numpy.abs(left.imag - right.imag) / imag_scale)

    finite_values = numpy.concatenate((real_mismatch[numpy.isfinite(real_mismatch)], imag_mismatch[numpy.isfinite(imag_mismatch)],))

    maximum = float(numpy.max(finite_values)) if finite_values.size else numpy.nan
    return (
        orders,
        real_left,
        real_right,
        imag_left,
        imag_right,
        real_mismatch,
        imag_mismatch,
        maximum,
    )


def _reduced_root_set_mismatch(reference: numpy.ndarray, candidate: numpy.ndarray) -> float:
    maximum = 0.0
    for reference_row, candidate_row in zip(reference, candidate):
        distances = numpy.abs(reference_row[:, None] - candidate_row[None, :])
        row_mismatch = max(float(numpy.max(numpy.min(distances, axis=1))), float(numpy.max(numpy.min(distances, axis=0))),)
        maximum = max(maximum, row_mismatch)
    return maximum


def _reduced_magnitude_order_departure_count(
    roots: numpy.ndarray, scale: float
) -> tuple[int, int]:
    """Return diagnostic departures from decreasing |Re(omega)| order.

    The counts do not modify the reduced roots.  A pointwise exchange would be
    admissible only at a true coalescence; otherwise it would replace two
    continuous analytic sheets by discontinuous index curves.
    """

    tolerance = _CANONICAL_REAL_ORDER_TOL * scale
    g2 = numpy.abs(roots[:, 2].real) + tolerance < numpy.abs(roots[:, 3].real)
    g3 = numpy.abs(roots[:, 4].real) + tolerance < numpy.abs(roots[:, 5].real)
    return int(numpy.count_nonzero(g2)), int(numpy.count_nonzero(g3))


def solve_six_reduced_branches(
    k: numpy.ndarray,
    wpe_negative: float,
    vd_negative: float,
    wpe_positive: float,
    vd_positive: float,
    wpi: float,
    *,
    slope_order_count: int = 4,
) -> SixComplexRoots:
    """Solve six branches on matched signed half-axes.

    The grid must be strictly increasing, contain one exact zero, and have
    matched positive/negative magnitudes. Both open half-axes are solved from
    their own polynomial coefficients and continued from small to large |k|.
    """

    _validate_species_order(vd_negative, vd_positive)
    k = numpy.asarray(k, dtype=float)
    if k.ndim != 1 or k.size < 3:
        raise ValueError("k must be one-dimensional and contain both signs.")
    if numpy.any(numpy.diff(k) <= 0.0):
        raise ValueError("k must be strictly increasing.")

    zero_indices = numpy.flatnonzero(k == 0.0)
    if zero_indices.size != 1:
        raise ValueError("k must contain exactly one exact zero.")
    zero_index = int(zero_indices[0])

    negative_indices = numpy.flatnonzero(k < 0.0)
    positive_indices = numpy.flatnonzero(k > 0.0)
    if negative_indices.size == 0 or positive_indices.size == 0:
        raise ValueError("k must contain nonzero values of both signs.")
    if negative_indices.size != positive_indices.size:
        raise ValueError("The signed half-axes must have equal lengths.")
    if not numpy.allclose(
        -k[negative_indices][::-1],
        k[positive_indices],
        rtol=5.0e-13,
        atol=0.0,
    ):
        raise ValueError("The grid must be matched under k -> -k.")

    # Direct independent continuations. The negative path is explicitly
    # reordered as -dk,-2dk,... so it also proceeds outward from k=0.
    positive_path = k[positive_indices]
    negative_path = k[negative_indices][::-1]
    positive = solve_reduced_half_axis(
        positive_path,
        wpe_negative,
        vd_negative,
        wpe_positive,
        vd_positive,
        wpi,
        sign=+1,
        slope_order_count=slope_order_count,
    )
    negative = solve_reduced_half_axis(
        negative_path,
        wpe_negative,
        vd_negative,
        wpe_positive,
        vd_positive,
        wpi,
        sign=-1,
        slope_order_count=slope_order_count,
    )

    roots = numpy.full((k.size, 6), numpy.nan + 1j * numpy.nan)
    raw_roots = numpy.full_like(roots, numpy.nan + 1j * numpy.nan)
    raw_source_index = numpy.full((k.size, 6), -1, dtype=int)
    solved_directly = numpy.zeros((k.size, 6), dtype=bool)
    conjugate_enforced = numpy.zeros(k.size, dtype=bool)
    tie_used = numpy.zeros(k.size, dtype=bool)

    roots[positive_indices] = positive.roots
    roots[negative_indices] = negative.roots[::-1]
    roots[zero_index] = reduced_limiting_roots_at_zero(wpe_negative, wpe_positive, wpi)

    raw_roots[positive_indices] = positive.raw_roots
    raw_roots[negative_indices] = negative.raw_roots[::-1]
    raw_source_index[positive_indices] = positive.raw_source_index
    raw_source_index[negative_indices] = negative.raw_source_index[::-1]
    solved_directly[positive_indices] = True
    solved_directly[negative_indices] = True
    conjugate_enforced[positive_indices] = positive.conjugate_partition_enforced
    conjugate_enforced[negative_indices] = (negative.conjugate_partition_enforced[::-1])
    tie_used[positive_indices] = positive.analytic_tie_break_used
    tie_used[negative_indices] = negative.analytic_tie_break_used[::-1]

    relative_residual = numpy.full((k.size, 6), numpy.nan, dtype=float)
    for index in numpy.concatenate((negative_indices, positive_indices)):
        for branch in range(6):
            relative_residual[index, branch] = reduced_relative_residual(
                k[index],
                roots[index, branch],
                wpe_negative,
                vd_negative,
                wpe_positive,
                vd_positive,
                wpi,
            )

    scale = _frequency_scale(wpe_negative, wpe_positive, wpi)
    same_real, opposite_imag, complex_pair_mask = _pair_consistency(roots, scale)

    # Separate branch-point estimates for G2 and G3 on each signed half-axis.
    k_b1_negative = _contiguous_complex_boundary(k[negative_indices], complex_pair_mask[negative_indices, 0])
    k_b1_positive = _contiguous_complex_boundary(k[positive_indices], complex_pair_mask[positive_indices, 0])
    k_b2_negative = _contiguous_complex_boundary(k[negative_indices], complex_pair_mask[negative_indices, 1])
    k_b2_positive = _contiguous_complex_boundary(k[positive_indices], complex_pair_mask[positive_indices, 1])
    def matched_boundary(first: float, second: float) -> float:
        if numpy.isfinite(first) and numpy.isfinite(second):
            return min(first, second)
        return numpy.nan

    k_b1_common = matched_boundary(k_b1_negative, k_b1_positive)
    k_b2_common = matched_boundary(k_b2_negative, k_b2_positive)

    slope_boundary_k = numpy.asarray([-k_b1_negative, +k_b1_positive, -k_b2_negative, +k_b2_positive], dtype=float,)
    (
        slope_orders,
        slope_real_left,
        slope_real_right,
        slope_imag_left,
        slope_imag_right,
        slope_real_relative_mismatch,
        slope_imag_relative_mismatch,
        maximum_rule5_relative_mismatch,
    ) = _reduced_rule5_slope_diagnostics(k, roots, slope_boundary_k, slope_order_count)

    # Validate the requested opposite-k identities without using them in the
    # signed-k solve. Rows are aligned by increasing positive |k|.
    aligned_negative = roots[negative_indices][::-1]
    opposite_target = -numpy.conjugate(roots[positive_indices][:, OPPOSITE_K_PARTNER])
    opposite_k_conversion_error = numpy.max(numpy.abs(aligned_negative - opposite_target), axis=0)
    opposite_k_root_set_mismatch = _reduced_root_set_mismatch(opposite_target, aligned_negative)

    selected_raw = numpy.take_along_axis(raw_roots[k != 0.0], raw_source_index[k != 0.0], axis=1)
    complete_root_assignment_error = float(numpy.max(numpy.abs(roots[k != 0.0] - selected_raw)))

    # Compare group growth magnitudes only where both pairs remain complex.
    nonzero = k != 0.0
    finite_common_boundaries = [value for value in (k_b1_common, k_b2_common) if numpy.isfinite(value)]
    common_limit = min(finite_common_boundaries) if finite_common_boundaries else numpy.inf
    growth_hierarchy_mask = (nonzero & complex_pair_mask[:, 0] & complex_pair_mask[:, 1] & (numpy.abs(k) < common_limit))
    g2_growth = 0.5 * (numpy.abs(roots[:, 2].imag) + numpy.abs(roots[:, 3].imag))
    g3_growth = 0.5 * (numpy.abs(roots[:, 4].imag) + numpy.abs(roots[:, 5].imag))
    growth_hierarchy_margin = (g2_growth - g3_growth) / scale
    if numpy.any(growth_hierarchy_mask):
        minimum_growth_hierarchy_margin = float(numpy.min(growth_hierarchy_margin[growth_hierarchy_mask]))
        growth_hierarchy_violation_count = int(numpy.count_nonzero(growth_hierarchy_margin[growth_hierarchy_mask] <= _GROWTH_COMPARISON_TOL))
    else:
        minimum_growth_hierarchy_margin = numpy.nan
        growth_hierarchy_violation_count = 0

    # Pointwise magnitude order is diagnostic only; analytic continuity has
    # priority over an index exchange away from coalescence.
    _reduced_magnitude_order_departure_count(roots, scale)

    k_bs = numpy.asarray(
        [
            [k_b1_negative, k_b1_positive, k_b1_common],
            [k_b2_negative, k_b2_positive, k_b2_common],
        ],
        dtype=float,
    )
    parameters = {
        "model": "reduced_cold_three_species",
        "wpe_negative": float(wpe_negative),
        "vd_negative": float(vd_negative),
        "wpe_positive": float(wpe_positive),
        "vd_positive": float(vd_positive),
        "wpi": float(wpi),
        "slope_order_count": int(slope_order_count),
    }

    return SixComplexRoots(
        k=k.copy(),
        roots=roots,
        k_bs=k_bs,
        parameters=parameters,
        raw_roots=raw_roots,
        raw_source_index=raw_source_index,
        solved_directly=solved_directly,
        relative_residual=relative_residual,
        same_real_error=same_real,
        opposite_imag_error=opposite_imag,
        complex_pair_mask=complex_pair_mask,
        conjugate_partition_enforced=conjugate_enforced,
        analytic_tie_break_used=tie_used,
        opposite_k_conversion_error=opposite_k_conversion_error,
        opposite_k_root_set_mismatch=float(opposite_k_root_set_mismatch),
        complete_root_assignment_error=complete_root_assignment_error,
        k_b1_negative=k_b1_negative,
        k_b1_positive=k_b1_positive,
        k_b1_common=k_b1_common,
        k_b2_negative=k_b2_negative,
        k_b2_positive=k_b2_positive,
        k_b2_common=k_b2_common,
        growth_hierarchy_mask=growth_hierarchy_mask,
        growth_hierarchy_margin=growth_hierarchy_margin,
        minimum_growth_hierarchy_margin=minimum_growth_hierarchy_margin,
        growth_hierarchy_violation_count=growth_hierarchy_violation_count,
        slope_boundary_k=slope_boundary_k,
        slope_orders=slope_orders,
        slope_real_left=slope_real_left,
        slope_real_right=slope_real_right,
        slope_imag_left=slope_imag_left,
        slope_imag_right=slope_imag_right,
        slope_real_relative_mismatch=slope_real_relative_mismatch,
        slope_imag_relative_mismatch=slope_imag_relative_mismatch,
        maximum_rule5_relative_mismatch=maximum_rule5_relative_mismatch,
    )

############################################################
# Full-Maxwellian continuation
############################################################
SQRT2 = float(numpy.sqrt(2.0))
SQRT_PI = float(numpy.sqrt(numpy.pi))

# Coefficients in 1 + eta Z(eta) ~ -sum_n a_n eta^(-2n).
_Q_COEFF = numpy.asarray(
    [
        1.0 / 2.0,
        3.0 / 4.0,
        15.0 / 8.0,
        105.0 / 16.0,
        945.0 / 32.0,
        10395.0 / 64.0,
        135135.0 / 128.0,
        2027025.0 / 256.0,
    ],
    dtype=float,
)


@dataclass(frozen=True)
class MaxwellianParameters:
    """Dimensionless cold-ion/two-electron parameters normalized by omega_pe and c."""

    alpha: numpy.ndarray
    drift_over_c: numpy.ndarray
    thermal_over_c: numpy.ndarray
    species_names: tuple[str, str, str]
    omega_scale: float


@dataclass(frozen=True)
class LocalRoot:
    """One corrected normalized root W=omega/omega_scale."""

    value: complex
    relative_residual: float
    converged: bool
    iterations: int


@dataclass(frozen=True)
class ContourResult:
    """Argument-principle result for a circular contour."""

    count: complex
    moments: numpy.ndarray
    minimum_boundary_value: float
    points: int
    valid: bool


@dataclass(frozen=True)
class HomotopyResult:
    """Six roots connected from lambda=0 to lambda=1 at one anchor."""

    roots: numpy.ndarray
    residual: numpy.ndarray
    converged: numpy.ndarray
    lambda_history: numpy.ndarray
    root_history: numpy.ndarray
    pseudoarclength_uses: int
    rejected_steps: int
    contour_certified: numpy.ndarray
    contour_count: numpy.ndarray


@dataclass(frozen=True)
class HalfAxisFullSolution:
    """Direct six-root solution on one open signed-k half-axis."""

    k_path: numpy.ndarray
    roots: numpy.ndarray
    reduced_seed: numpy.ndarray
    relative_residual: numpy.ndarray
    converged: numpy.ndarray
    iterations: numpy.ndarray
    candidate_count: numpy.ndarray
    pseudoarclength_used: numpy.ndarray
    contour_used: numpy.ndarray
    contour_certified: numpy.ndarray
    contour_count: numpy.ndarray
    analytic_tie_break_used: numpy.ndarray
    anchor_index: int
    homotopy: HomotopyResult




############################################################
# Fried--Conte response and dimensionless dispersion functions
############################################################


def _response_and_derivative(eta: complex) -> tuple[complex, complex]:
    """Return q=1+eta Z(eta) and dq/deta on the causal Landau sheet."""

    eta = complex(eta)
    if eta == 0.0:
        return 1.0 + 0.0j, 1j * SQRT_PI

    eta2 = eta * eta
    # In this sector the exponentially small Landau term is negligible and
    # the algebraic expansion avoids cancellation in 1+eta*Z.
    if abs(eta) >= 18.0 and eta2.real >= 40.0:
        inverse_square = 1.0 / eta2
        power = inverse_square
        response = 0.0 + 0.0j
        derivative = 0.0 + 0.0j
        for order, coefficient in enumerate(_Q_COEFF, start=1):
            response -= coefficient * power
            derivative += 2.0 * order * coefficient * power / eta
            power *= inverse_square
        return complex(response), complex(derivative)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with numpy.errstate(over="ignore", invalid="ignore", under="ignore"):
            z_value = complex(1j * SQRT_PI * wofz(eta))
            response = 1.0 + eta * z_value
            derivative = z_value - 2.0 * eta * response
    return complex(response), complex(derivative)


def build_maxwellian_parameters(
    wpe_negative: float,
    vd_negative: float,
    vthe_negative: float,
    wpe_positive: float,
    vd_positive: float,
    vthe_positive: float,
    wpi: float,
    vthi: float,
) -> MaxwellianParameters:
    """Construct normalized ion/electron parameters."""

    values = numpy.asarray(
        [
            wpe_negative,
            vd_negative,
            vthe_negative,
            wpe_positive,
            vd_positive,
            vthe_positive,
            wpi,
            vthi,
        ],
        dtype=float,
    )
    if not numpy.isfinite(values).all():
        raise ValueError("All plasma parameters must be finite.")
    if min(wpe_negative, wpe_positive, wpi) <= 0.0:
        raise ValueError("All plasma frequencies must be positive.")
    if min(vthe_negative, vthe_positive, vthi) <= 0.0:
        raise ValueError("All thermal speeds must be positive.")
    if not vd_negative < 0.0 or not vd_positive > 0.0:
        raise ValueError("Require vd_negative<0 and vd_positive>0.")

    omega_scale = float(numpy.hypot(wpe_negative, wpe_positive))
    return MaxwellianParameters(
        alpha=numpy.asarray(
            [
                (wpi / omega_scale) ** 2,
                (wpe_negative / omega_scale) ** 2,
                (wpe_positive / omega_scale) ** 2,
            ],
            dtype=float,
        ),
        drift_over_c=numpy.asarray(
            [0.0, vd_negative / CGS["c"], vd_positive / CGS["c"]],
            dtype=float,
        ),
        thermal_over_c=numpy.asarray(
            [vthi / CGS["c"], vthe_negative / CGS["c"], vthe_positive / CGS["c"]],
            dtype=float,
        ),
        species_names=("ion", "electron_negative", "electron_positive"),
        omega_scale=omega_scale,
    )


def maxwellian_scaled_equation(
    K: float,
    W: complex,
    parameters: MaxwellianParameters,
) -> tuple[complex, complex, complex, float]:
    """Return F=K^2 D for cold ions and two Maxwellian electrons.

    The stationary ion population is retained through its cold inertial term,
    while both electron populations use the causal Maxwellian response.  This
    model has the exact sixth-order reduced equation as its zero-temperature
    limit and avoids introducing unrelated ion Landau sheets into the selected
    electron/electron--ion branch continuation.
    """

    if K == 0.0:
        nan = numpy.nan + 1j * numpy.nan
        return nan, nan, nan, numpy.nan

    sign = 1.0 if K > 0.0 else -1.0
    absolute_K = abs(K)
    alpha_i = float(parameters.alpha[0])
    if abs(W) <= 1.0e-300:
        nan = numpy.nan + 1j * numpy.nan
        return nan, nan, nan, numpy.inf

    ion_term = -alpha_i * K * K / (W * W)
    value = complex(K * K + ion_term)
    derivative_W = 2.0 * alpha_i * K * K / (W**3)
    derivative_K = 2.0 * K - 2.0 * alpha_i * K / (W * W)
    scale = float(K * K + abs(ion_term))

    with numpy.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for alpha, drift, thermal in zip(
            parameters.alpha[1:],
            parameters.drift_over_c[1:],
            parameters.thermal_over_c[1:],
        ):
            eta = (W - K * drift) / (SQRT2 * absolute_K * thermal)
            response, response_derivative = _response_and_derivative(eta)
            coefficient = alpha / (thermal * thermal)
            term = coefficient * response
            value += term
            scale += abs(term)

            derivative_eta_W = 1.0 / (SQRT2 * absolute_K * thermal)
            derivative_eta_K = -sign * W / (SQRT2 * thermal * K * K)
            derivative_W += coefficient * response_derivative * derivative_eta_W
            derivative_K += coefficient * response_derivative * derivative_eta_K

    return complex(value), complex(derivative_W), complex(derivative_K), float(scale)

def cold_scaled_equation(
    K: float,
    W: complex,
    parameters: MaxwellianParameters,
) -> tuple[complex, complex, complex, float]:
    """Return the K^2-scaled sixth-order cold rational equation."""

    if K == 0.0:
        nan = numpy.nan + 1j * numpy.nan
        return nan, nan, nan, numpy.nan

    value = complex(K * K)
    derivative_W = 0.0 + 0.0j
    derivative_K = complex(2.0 * K)
    scale = float(K * K)

    for alpha, drift in zip(parameters.alpha, parameters.drift_over_c):
        delta = W - K * drift
        if abs(delta) <= 1.0e-300:
            nan = numpy.nan + 1j * numpy.nan
            return nan, nan, nan, numpy.inf
        term = -alpha * K * K / (delta * delta)
        value += term
        scale += abs(term)
        derivative_W += 2.0 * alpha * K * K / (delta**3)
        derivative_K += -alpha * (
            2.0 * K / (delta * delta)
            + 2.0 * K * K * drift / (delta**3)
        )

    return complex(value), complex(derivative_W), complex(derivative_K), float(scale)


def homotopy_scaled_equation(
    K: float,
    W: complex,
    lam: float,
    parameters: MaxwellianParameters,
) -> tuple[complex, complex, complex, complex, float]:
    """Thermal homotopy for the two electron populations.

    The ion inertial term remains cold for all lambda.  Electron thermal speeds
    are scaled as v_th,e(lambda)=sqrt(lambda)v_th,e.  The lambda->0 limit is
    therefore exactly the uploaded sixth-order reduced equation.
    """

    if lam < 0.0:
        nan = numpy.nan + 1j * numpy.nan
        return nan, nan, nan, nan, numpy.nan
    if lam <= 1.0e-14:
        cold_value, cold_W, cold_K, cold_scale = cold_scaled_equation(
            K, W, parameters
        )
        derivative_lam = 0.0 + 0.0j
        for alpha, drift, thermal in zip(
            parameters.alpha[1:],
            parameters.drift_over_c[1:],
            parameters.thermal_over_c[1:],
        ):
            delta = W - K * drift
            if abs(delta) <= 1.0e-300:
                nan = numpy.nan + 1j * numpy.nan
                return nan, nan, nan, nan, numpy.inf
            derivative_lam -= 3.0 * alpha * K**4 * thermal**2 / (delta**4)
        return cold_value, cold_W, cold_K, complex(derivative_lam), cold_scale

    sign = 1.0 if K > 0.0 else -1.0
    absolute_K = abs(K)
    square_root_lam = float(numpy.sqrt(lam))
    alpha_i = float(parameters.alpha[0])
    if abs(W) <= 1.0e-300:
        nan = numpy.nan + 1j * numpy.nan
        return nan, nan, nan, nan, numpy.inf

    ion_term = -alpha_i * K * K / (W * W)
    value = complex(K * K + ion_term)
    derivative_W = 2.0 * alpha_i * K * K / (W**3)
    derivative_K = 2.0 * K - 2.0 * alpha_i * K / (W * W)
    derivative_lam = 0.0 + 0.0j
    scale = float(K * K + abs(ion_term))

    with numpy.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for alpha, drift, thermal in zip(
            parameters.alpha[1:],
            parameters.drift_over_c[1:],
            parameters.thermal_over_c[1:],
        ):
            effective_thermal = square_root_lam * thermal
            eta = (W - K * drift) / (
                SQRT2 * absolute_K * effective_thermal
            )
            response, response_derivative = _response_and_derivative(eta)
            coefficient = alpha / (lam * thermal * thermal)
            term = coefficient * response
            value += term
            scale += abs(term)

            derivative_eta_W = 1.0 / (
                SQRT2 * absolute_K * effective_thermal
            )
            derivative_eta_K = -sign * W / (
                SQRT2 * effective_thermal * K * K
            )
            derivative_W += coefficient * response_derivative * derivative_eta_W
            derivative_K += coefficient * response_derivative * derivative_eta_K
            derivative_lam -= (
                alpha
                / (thermal * thermal * lam * lam)
                * (response + 0.5 * eta * response_derivative)
            )

    return (
        complex(value),
        complex(derivative_W),
        complex(derivative_K),
        complex(derivative_lam),
        float(scale),
    )

def full_maxwellian_relative_residual(
    k: float,
    omega: complex,
    parameters: MaxwellianParameters,
) -> float:
    if k == 0.0:
        return numpy.nan
    K = k * CGS["c"] / parameters.omega_scale
    W = omega / parameters.omega_scale
    value, _, _, scale = maxwellian_scaled_equation(K, W, parameters)
    if not _finite_complex(value) or not numpy.isfinite(scale) or scale <= 0.0:
        return numpy.inf
    return float(abs(value) / scale)


############################################################
# Local Newton and pseudo-arclength correctors
############################################################


def _finite_complex(value: complex) -> bool:
    return bool(numpy.isfinite(value.real) and numpy.isfinite(value.imag))


def _complex_jacobian(derivative: complex) -> numpy.ndarray:
    return numpy.asarray(
        [
            [derivative.real, -derivative.imag],
            [derivative.imag, derivative.real],
        ],
        dtype=float,
    )


def _solve_fixed_complex_root(
    seed: complex,
    evaluator: Callable[[complex], tuple[complex, complex, float]],
    *,
    residual_tolerance: float,
    maximum_iterations: int = 45,
    maximum_step: float = 0.35,
) -> LocalRoot:
    """Correct one complex seed with analytic Newton and a real fallback."""

    current = complex(seed)
    best = current
    best_residual = numpy.inf

    for iteration in range(1, maximum_iterations + 1):
        value, derivative, scale = evaluator(current)
        if (
            not _finite_complex(value)
            or not _finite_complex(derivative)
            or not numpy.isfinite(scale)
            or scale <= 0.0
            or abs(derivative) <= 1.0e-16
        ):
            break
        residual = float(abs(value) / scale)
        if residual < best_residual:
            best = current
            best_residual = residual
        if residual <= residual_tolerance:
            return LocalRoot(current, residual, True, iteration)

        step = value / derivative
        local_limit = maximum_step * max(1.0, abs(current))
        if abs(step) > local_limit:
            step *= local_limit / abs(step)

        accepted = False
        for backtrack in range(14):
            trial = current - step * (0.5**backtrack)
            trial_value, _, trial_scale = evaluator(trial)
            if (
                _finite_complex(trial_value)
                and numpy.isfinite(trial_scale)
                and trial_scale > 0.0
                and abs(trial_value) / trial_scale < residual
            ):
                current = trial
                accepted = True
                break
        if not accepted:
            current -= 0.03 * step

        if abs(step) <= 2.0e-13 * max(1.0, abs(current)):
            break

    def real_system(vector: numpy.ndarray) -> numpy.ndarray:
        trial = complex(vector[0], vector[1])
        value, _, scale = evaluator(trial)
        if not _finite_complex(value) or not numpy.isfinite(scale) or scale <= 0.0:
            return numpy.asarray([1.0e100, 1.0e100], dtype=float)
        return numpy.asarray([value.real / scale, value.imag / scale], dtype=float)

    def real_jacobian(vector: numpy.ndarray) -> numpy.ndarray:
        trial = complex(vector[0], vector[1])
        _, derivative, scale = evaluator(trial)
        if not _finite_complex(derivative) or not numpy.isfinite(scale) or scale <= 0.0:
            return numpy.eye(2) * 1.0e100
        return _complex_jacobian(derivative) / scale

    try:
        fallback = root(
            real_system,
            numpy.asarray([best.real, best.imag], dtype=float),
            jac=real_jacobian,
            method="hybr",
            options={"xtol": 1.0e-11, "maxfev": 350},
        )
        trial = complex(fallback.x[0], fallback.x[1])
        value, _, scale = evaluator(trial)
        trial_residual = (
            float(abs(value) / scale)
            if _finite_complex(value) and numpy.isfinite(scale) and scale > 0.0
            else numpy.inf
        )
        if trial_residual < best_residual:
            best = trial
            best_residual = trial_residual
    except Exception:
        pass

    return LocalRoot(
        best,
        float(best_residual),
        bool(best_residual <= residual_tolerance),
        maximum_iterations,
    )


def solve_full_root_at_K(
    K: float,
    seed_W: complex,
    parameters: MaxwellianParameters,
    *,
    residual_tolerance: float = 5.0e-11,
) -> LocalRoot:
    if K == 0.0:
        raise ValueError("The Maxwellian equation is not evaluated at K=0.")

    def evaluator(W: complex) -> tuple[complex, complex, float]:
        value, derivative_W, _, scale = maxwellian_scaled_equation(
            K, W, parameters
        )
        return value, derivative_W, scale

    return _solve_fixed_complex_root(
        seed_W,
        evaluator,
        residual_tolerance=residual_tolerance,
    )


def solve_homotopy_root_at_lambda(
    K: float,
    lam: float,
    seed_W: complex,
    parameters: MaxwellianParameters,
    *,
    residual_tolerance: float,
) -> LocalRoot:
    def evaluator(W: complex) -> tuple[complex, complex, float]:
        value, derivative_W, _, _, scale = homotopy_scaled_equation(
            K, W, lam, parameters
        )
        return value, derivative_W, scale

    return _solve_fixed_complex_root(
        seed_W,
        evaluator,
        residual_tolerance=residual_tolerance,
    )


def _normalized_tangent(vector: numpy.ndarray, preferred_last_sign: float) -> numpy.ndarray:
    norm = float(numpy.linalg.norm(vector))
    if not numpy.isfinite(norm) or norm <= 1.0e-16:
        result = numpy.zeros_like(vector, dtype=float)
        result[-1] = preferred_last_sign
        return result
    result = numpy.asarray(vector, dtype=float) / norm
    if result[-1] * preferred_last_sign < 0.0:
        result = -result
    return result


def _pseudoarclength_homotopy_seed(
    K: float,
    current_lam: float,
    current_W: complex,
    previous_lam: float | None,
    previous_W: complex | None,
    target_lam: float,
    parameters: MaxwellianParameters,
) -> tuple[complex, bool]:
    """Return a pseudo-arclength-assisted seed at a target lambda."""

    if previous_lam is not None and previous_W is not None:
        tangent = _normalized_tangent(
            numpy.asarray(
                [
                    current_W.real - previous_W.real,
                    current_W.imag - previous_W.imag,
                    current_lam - previous_lam,
                ],
                dtype=float,
            ),
            +1.0,
        )
    else:
        value, derivative_W, _, derivative_lam, _ = homotopy_scaled_equation(
            K, current_W, current_lam, parameters
        )
        del value
        if abs(derivative_W) <= 1.0e-16 or not _finite_complex(derivative_lam):
            derivative_path = 0.0 + 0.0j
        else:
            derivative_path = -derivative_lam / derivative_W
        tangent = _normalized_tangent(
            numpy.asarray(
                [derivative_path.real, derivative_path.imag, 1.0], dtype=float
            ),
            +1.0,
        )

    if abs(tangent[2]) <= 1.0e-8:
        predicted = numpy.asarray(
            [current_W.real, current_W.imag, target_lam], dtype=float
        )
    else:
        arc_step = (target_lam - current_lam) / tangent[2]
        predicted = numpy.asarray(
            [current_W.real, current_W.imag, current_lam], dtype=float
        ) + arc_step * tangent

    def system(vector: numpy.ndarray) -> numpy.ndarray:
        W = complex(vector[0], vector[1])
        lam = float(vector[2])
        value, _, _, _, scale = homotopy_scaled_equation(K, W, lam, parameters)
        if not _finite_complex(value) or not numpy.isfinite(scale) or scale <= 0.0:
            return numpy.asarray([1.0e100, 1.0e100, 1.0e100])
        return numpy.asarray(
            [
                value.real / scale,
                value.imag / scale,
                float(numpy.dot(tangent, vector - predicted)),
            ],
            dtype=float,
        )

    def jacobian(vector: numpy.ndarray) -> numpy.ndarray:
        W = complex(vector[0], vector[1])
        lam = float(vector[2])
        _, derivative_W, _, derivative_lam, scale = homotopy_scaled_equation(
            K, W, lam, parameters
        )
        matrix = numpy.zeros((3, 3), dtype=float)
        matrix[0:2, 0:2] = _complex_jacobian(derivative_W) / scale
        matrix[0, 2] = derivative_lam.real / scale
        matrix[1, 2] = derivative_lam.imag / scale
        matrix[2] = tangent
        return matrix

    try:
        answer = root(
            system,
            predicted,
            jac=jacobian,
            method="hybr",
            options={"xtol": 2.0e-11, "maxfev": 250},
        )
        seed = complex(answer.x[0], answer.x[1])
        valid = bool(
            answer.success
            and _finite_complex(seed)
            and -0.15 <= answer.x[2] <= 1.15
        )
        return seed, valid
    except Exception:
        return complex(predicted[0], predicted[1]), False


def _pseudoarclength_K_seed(
    previous_K: float,
    previous_W: complex,
    earlier_K: float,
    earlier_W: complex,
    target_K: float,
    parameters: MaxwellianParameters,
) -> tuple[complex, bool]:
    """Return a pseudo-arclength-assisted seed for a new signed K."""

    preferred = 1.0 if target_K > previous_K else -1.0
    tangent = _normalized_tangent(
        numpy.asarray(
            [
                previous_W.real - earlier_W.real,
                previous_W.imag - earlier_W.imag,
                previous_K - earlier_K,
            ],
            dtype=float,
        ),
        preferred,
    )
    if abs(tangent[2]) <= 1.0e-8:
        predicted = numpy.asarray(
            [previous_W.real, previous_W.imag, target_K], dtype=float
        )
    else:
        arc_step = (target_K - previous_K) / tangent[2]
        predicted = numpy.asarray(
            [previous_W.real, previous_W.imag, previous_K], dtype=float
        ) + arc_step * tangent

    def system(vector: numpy.ndarray) -> numpy.ndarray:
        W = complex(vector[0], vector[1])
        K = float(vector[2])
        if K == 0.0 or K * target_K <= 0.0:
            return numpy.asarray([1.0e100, 1.0e100, 1.0e100])
        value, _, _, scale = maxwellian_scaled_equation(K, W, parameters)
        if not _finite_complex(value) or not numpy.isfinite(scale) or scale <= 0.0:
            return numpy.asarray([1.0e100, 1.0e100, 1.0e100])
        return numpy.asarray(
            [
                value.real / scale,
                value.imag / scale,
                float(numpy.dot(tangent, vector - predicted)),
            ],
            dtype=float,
        )

    def jacobian(vector: numpy.ndarray) -> numpy.ndarray:
        W = complex(vector[0], vector[1])
        K = float(vector[2])
        _, derivative_W, derivative_K, scale = maxwellian_scaled_equation(
            K, W, parameters
        )
        matrix = numpy.zeros((3, 3), dtype=float)
        matrix[0:2, 0:2] = _complex_jacobian(derivative_W) / scale
        matrix[0, 2] = derivative_K.real / scale
        matrix[1, 2] = derivative_K.imag / scale
        matrix[2] = tangent
        return matrix

    try:
        answer = root(
            system,
            predicted,
            jac=jacobian,
            method="hybr",
            options={"xtol": 2.0e-11, "maxfev": 220},
        )
        seed = complex(answer.x[0], answer.x[1])
        valid = bool(
            answer.success
            and _finite_complex(seed)
            and answer.x[2] * target_K > 0.0
        )
        return seed, valid
    except Exception:
        return complex(predicted[0], predicted[1]), False


############################################################
# Candidate management and strict complete-root group assignment
############################################################


def _append_unique_candidate(
    candidates: list[LocalRoot],
    candidate: LocalRoot,
    separation: float,
) -> None:
    if not candidate.converged or not _finite_complex(candidate.value):
        return
    for index, existing in enumerate(candidates):
        if abs(candidate.value - existing.value) <= separation:
            if candidate.relative_residual < existing.relative_residual:
                candidates[index] = candidate
            return
    candidates.append(candidate)


def _pair_descriptor(pair: numpy.ndarray) -> tuple[complex, complex]:
    center = 0.5 * (pair[0] + pair[1])
    half = 0.5 * (pair[0] - pair[1])
    return center, half * half


def _pair_set_distance(first: numpy.ndarray, second: numpy.ndarray) -> float:
    direct = abs(first[0] - second[0]) + abs(first[1] - second[1])
    crossed = abs(first[0] - second[1]) + abs(first[1] - second[0])
    return float(min(direct, crossed))


def _order_pair_by_continuity(
    pair: numpy.ndarray,
    residual: numpy.ndarray,
    iterations: numpy.ndarray,
    predicted: numpy.ndarray,
    reduced: numpy.ndarray,
    sign: int,
    group_index: int,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, bool]:
    """Order an intact pair by complete-root continuity.

    The direct and crossed assignments are compared in the complex-omega
    plane.  A deterministic convention is used only when these assignments are
    numerically degenerate, as at a coalescence or analytic reconnection.  No
    pointwise sorting by Re(omega), Im(omega), or |Re(omega)| is applied after a
    nondegenerate continuity decision.
    """

    direct = abs(pair[0] - predicted[0]) + abs(pair[1] - predicted[1])
    crossed = abs(pair[1] - predicted[0]) + abs(pair[0] - predicted[1])
    reduced_direct = abs(pair[0] - reduced[0]) + abs(pair[1] - reduced[1])
    reduced_crossed = abs(pair[1] - reduced[0]) + abs(pair[0] - reduced[1])
    direct += 0.02 * reduced_direct
    crossed += 0.02 * reduced_crossed

    scale = max(1.0, abs(pair[0]), abs(pair[1]), abs(predicted[0]), abs(predicted[1]))
    tolerance = 2.0e-7 * scale + 2.0e-5 * min(direct, crossed)
    tie = abs(direct - crossed) <= tolerance

    if not tie:
        order = numpy.asarray([0, 1] if direct < crossed else [1, 0], dtype=int)
    else:
        # The convention below is invoked only where both complete-root
        # assignments are indistinguishable within numerical tolerance.
        if group_index in (1, 2) and numpy.max(numpy.abs(pair.imag)) > 2.0e-9:
            order = numpy.argsort(-pair.imag)
        elif group_index in (1, 2):
            order = numpy.argsort(sign * pair.real)
        else:
            order = numpy.argsort(-pair.real)

    return pair[order], residual[order], iterations[order], bool(tie)


def _assign_candidates_strict_groups(
    candidates: list[LocalRoot],
    predicted: numpy.ndarray,
    reduced: numpy.ndarray,
    sign: int,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Assign six candidates in the strict order G1, G2, G3."""

    if len(candidates) < 6:
        raise RuntimeError("Fewer than six distinct kinetic roots were found.")

    available = list(range(len(candidates)))
    assigned = numpy.empty(6, dtype=numpy.complex128)
    residual = numpy.empty(6, dtype=float)
    iterations = numpy.empty(6, dtype=int)
    ties = numpy.zeros(3, dtype=bool)

    for group_index, columns in enumerate(GROUP_COLUMNS):
        target_pred = predicted[list(columns)]
        target_reduced = reduced[list(columns)]
        best: tuple[float, int, int] | None = None

        pred_center, pred_q = _pair_descriptor(target_pred)
        reduced_center, reduced_q = _pair_descriptor(target_reduced)
        for first, second in combinations(available, 2):
            pair = numpy.asarray(
                [candidates[first].value, candidates[second].value],
                dtype=numpy.complex128,
            )
            center, q = _pair_descriptor(pair)
            cost = _pair_set_distance(pair, target_pred)
            # Branchwise complete-root continuity dominates.  Pair invariants
            # and reduced ancestry are only tie-breakers near coalescence.
            cost += 0.03 * _pair_set_distance(pair, target_reduced)
            cost += 0.03 * abs(center - pred_center)
            cost += 0.005 * abs(q - pred_q) / max(1.0, abs(pred_q))
            cost += 0.01 * abs(center - reduced_center)
            cost += 0.002 * abs(q - reduced_q) / max(1.0, abs(reduced_q))
            candidate_key = (float(cost), first, second)
            if best is None or candidate_key < best:
                best = candidate_key

        if best is None:
            raise RuntimeError(f"No admissible pair was found for group {group_index+1}.")

        indices = numpy.asarray([best[1], best[2]], dtype=int)
        pair = numpy.asarray([candidates[i].value for i in indices])
        pair_residual = numpy.asarray([candidates[i].relative_residual for i in indices])
        pair_iterations = numpy.asarray([candidates[i].iterations for i in indices])
        pair, pair_residual, pair_iterations, tie = _order_pair_by_continuity(
            pair,
            pair_residual,
            pair_iterations,
            target_pred,
            target_reduced,
            sign,
            group_index,
        )
        assigned[list(columns)] = pair
        residual[list(columns)] = pair_residual
        iterations[list(columns)] = pair_iterations
        ties[group_index] = tie
        for index in sorted(indices.tolist(), reverse=True):
            available.remove(index)

    # Near an inner four-root reconnection, the two pair centers coalesce and
    # the two squared-separation sheets split.  Pointwise nearest-neighbor
    # matching then follows the locally closest pair and can exchange the
    # physical G2/G3 ancestry on a fine k grid.  When both the predicted and
    # corrected pair centers are in the reconnection neighborhood, retain the
    # stronger-radius pair as G2 and the weaker-radius pair as G3.  The exchange
    # acts on the two complete complex pairs and is never applied away from a
    # center coalescence.
    first_inner = assigned[2:4].copy()
    second_inner = assigned[4:6].copy()
    first_center, _ = _pair_descriptor(first_inner)
    second_center, _ = _pair_descriptor(second_inner)
    predicted_first_center, _ = _pair_descriptor(predicted[2:4])
    predicted_second_center, _ = _pair_descriptor(predicted[4:6])
    first_radius = 0.5 * abs(first_inner[0] - first_inner[1])
    second_radius = 0.5 * abs(second_inner[0] - second_inner[1])
    predicted_first_radius = 0.5 * abs(predicted[2] - predicted[3])
    predicted_second_radius = 0.5 * abs(predicted[4] - predicted[5])
    current_radius_sum = max(first_radius + second_radius, 1.0e-10)
    predicted_radius_sum = max(
        predicted_first_radius + predicted_second_radius, 1.0e-10
    )
    current_centers_close = (
        abs(first_center - second_center) <= 0.35 * current_radius_sum
    )
    predicted_centers_close = (
        abs(predicted_first_center - predicted_second_center)
        <= 0.65 * predicted_radius_sum
    )
    radius_split = abs(first_radius - second_radius) > (
        2.0e-6 + 2.0e-4 * max(first_radius, second_radius, 1.0)
    )
    if (
        current_centers_close
        and predicted_centers_close
        and radius_split
        and second_radius > first_radius
    ):
        assigned[2:4], assigned[4:6] = (
            assigned[4:6].copy(),
            assigned[2:4].copy(),
        )
        residual[2:4], residual[4:6] = (
            residual[4:6].copy(),
            residual[2:4].copy(),
        )
        iterations[2:4], iterations[4:6] = (
            iterations[4:6].copy(),
            iterations[2:4].copy(),
        )
        ties[1] = True
        ties[2] = True

    return assigned, residual, iterations, ties


def _multiscale_history_predict(
    values: numpy.ndarray,
    grid: numpy.ndarray,
    history: list[int],
    target_index: int,
    maximum_order: int,
) -> numpy.ndarray:
    current_index = history[-1]
    current = values[current_index]
    if len(history) < 2 or maximum_order <= 0:
        return current.copy()

    predictions = []
    weights = []
    for order in range(1, min(maximum_order, len(history) - 1) + 1):
        earlier_index = history[-1 - order]
        denominator = grid[current_index] - grid[earlier_index]
        if denominator == 0.0:
            continue
        slope = (current - values[earlier_index]) / denominator
        predictions.append(current + slope * (grid[target_index] - grid[current_index]))
        weights.append(1.0 / order)
    if not predictions:
        return current.copy()
    weight_array = numpy.asarray(weights, dtype=float)
    return numpy.tensordot(
        weight_array / weight_array.sum(),
        numpy.stack(predictions, axis=0),
        axes=(0, 0),
    )


############################################################
# Argument-principle root counts and contour moments
############################################################


def argument_principle_circle(
    K: float,
    center: complex,
    radius: float,
    parameters: MaxwellianParameters,
    *,
    maximum_moment: int = 2,
    points: int = 128,
) -> ContourResult:
    """Count full-Maxwellian zeros and evaluate power sums inside a circle."""

    if K == 0.0 or radius <= 0.0 or points < 32:
        return ContourResult(
            numpy.nan + 1j * numpy.nan,
            numpy.full(maximum_moment + 1, numpy.nan + 1j * numpy.nan),
            numpy.nan,
            points,
            False,
        )

    theta = 2.0 * numpy.pi * numpy.arange(points, dtype=float) / points
    direction = numpy.exp(1j * theta)
    contour = center + radius * direction
    logarithmic_derivative = numpy.empty(points, dtype=numpy.complex128)
    minimum_value = numpy.inf

    for index, W in enumerate(contour):
        value, derivative_W, _, _ = maxwellian_scaled_equation(K, W, parameters)
        if not _finite_complex(value) or not _finite_complex(derivative_W) or value == 0.0:
            return ContourResult(
                numpy.nan + 1j * numpy.nan,
                numpy.full(maximum_moment + 1, numpy.nan + 1j * numpy.nan),
                0.0,
                points,
                False,
            )
        logarithmic_derivative[index] = derivative_W / value
        minimum_value = min(minimum_value, abs(value))

    # dW = i*r*exp(i theta)dtheta.  Division by 2*pi*i leaves the mean below.
    weight = radius * direction / points
    moments = numpy.empty(maximum_moment + 1, dtype=numpy.complex128)
    for order in range(maximum_moment + 1):
        moments[order] = numpy.sum(
            (contour**order) * logarithmic_derivative * weight
        )
    count = moments[0]
    nearest_integer = int(numpy.rint(count.real))
    valid = bool(
        numpy.isfinite(count.real)
        and numpy.isfinite(count.imag)
        and abs(count.imag) <= 0.08
        and abs(count.real - nearest_integer) <= 0.08
        and minimum_value > 1.0e-13
    )
    return ContourResult(
        complex(count), moments, float(minimum_value), int(points), valid
    )


def _adaptive_contour(
    K: float,
    center: complex,
    radius: float,
    expected_count: int,
    parameters: MaxwellianParameters,
    *,
    maximum_moment: int,
    base_points: int,
) -> ContourResult:
    best: ContourResult | None = None
    for factor in (1.0, 0.85, 1.15, 0.70, 1.30):
        trial_radius = radius * factor
        for points in (base_points, 2 * base_points, 4 * base_points):
            result = argument_principle_circle(
                K,
                center,
                trial_radius,
                parameters,
                maximum_moment=maximum_moment,
                points=points,
            )
            if best is None or (
                result.valid
                and abs(result.count.real - expected_count)
                < abs(best.count.real - expected_count)
            ):
                best = result
            if result.valid and abs(result.count.real - expected_count) <= 0.08:
                return result
    assert best is not None
    return best


def _contour_certify_and_recover(
    K: float,
    roots_in: numpy.ndarray,
    parameters: MaxwellianParameters,
    *,
    residual_tolerance: float,
    contour_points: int,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Certify each strict group and recover roots from contour moments."""

    roots_out = roots_in.copy()
    certified = numpy.zeros(6, dtype=bool)
    counts = numpy.full(6, numpy.nan, dtype=float)

    for columns in GROUP_COLUMNS:
        group = roots_out[list(columns)]
        outside = [index for index in range(6) if index not in columns]
        center, _ = _pair_descriptor(group)
        group_radius = max(abs(group[0] - center), abs(group[1] - center))
        outside_distance = min(abs(roots_out[index] - center) for index in outside)

        use_pair_contour = bool(
            group_radius > 0.0
            and 1.45 * group_radius < 0.72 * outside_distance
        )
        if use_pair_contour:
            radius = min(1.75 * group_radius + 1.0e-6, 0.70 * outside_distance)
            contour = _adaptive_contour(
                K,
                center,
                radius,
                2,
                parameters,
                maximum_moment=2,
                base_points=contour_points,
            )
            counts[list(columns)] = contour.count.real
            if contour.valid and abs(contour.count.real - 2.0) <= 0.08:
                p1 = contour.moments[1]
                p2 = contour.moments[2]
                recovered_center = 0.5 * p1
                recovered_q = 0.5 * p2 - recovered_center * recovered_center
                half = numpy.sqrt(recovered_q)
                recovered = numpy.asarray(
                    [recovered_center + half, recovered_center - half],
                    dtype=numpy.complex128,
                )
                polished: list[LocalRoot] = []
                for seed in recovered:
                    polished.append(
                        solve_full_root_at_K(
                            K,
                            seed,
                            parameters,
                            residual_tolerance=residual_tolerance,
                        )
                    )
                if all(item.converged for item in polished):
                    # Preserve the incoming complete-root order.
                    direct = abs(polished[0].value - group[0]) + abs(
                        polished[1].value - group[1]
                    )
                    crossed = abs(polished[1].value - group[0]) + abs(
                        polished[0].value - group[1]
                    )
                    order = (0, 1) if direct <= crossed else (1, 0)
                    roots_out[list(columns)] = numpy.asarray(
                        [polished[index].value for index in order]
                    )
                    certified[list(columns)] = True
                    continue

        # Fall back to one-root contours when the pair cannot be isolated.
        for branch in columns:
            nearest = min(
                abs(roots_out[branch] - roots_out[other])
                for other in range(6)
                if other != branch
            )
            radius = max(2.0e-5, 0.32 * nearest)
            contour = _adaptive_contour(
                K,
                roots_out[branch],
                radius,
                1,
                parameters,
                maximum_moment=1,
                base_points=contour_points,
            )
            counts[branch] = contour.count.real
            if contour.valid and abs(contour.count.real - 1.0) <= 0.08:
                seed = contour.moments[1]
                polished = solve_full_root_at_K(
                    K,
                    seed,
                    parameters,
                    residual_tolerance=residual_tolerance,
                )
                if polished.converged:
                    roots_out[branch] = polished.value
                    certified[branch] = True

    return roots_out, certified, counts


############################################################
# Anchor homotopy
############################################################


def _broad_seed_offsets(scale: float) -> tuple[complex, ...]:
    radial = max(0.004, 0.035 * scale)
    return (
        0.0j,
        +radial,
        -radial,
        +1j * radial,
        -1j * radial,
        +radial * (1.0 + 1.0j),
        +radial * (1.0 - 1.0j),
        -radial * (1.0 + 1.0j),
        -radial * (1.0 - 1.0j),
    )


def solve_anchor_homotopy(
    K: float,
    reduced_roots_W: numpy.ndarray,
    parameters: MaxwellianParameters,
    *,
    sign: int,
    residual_tolerance: float,
    initial_step: float = 0.04,
    minimum_step: float = 2.5e-4,
    maximum_step: float = 0.12,
    contour_points: int = 128,
) -> HomotopyResult:
    """Connect all six reduced roots to full-Maxwellian roots at one anchor."""

    reduced_roots_W = numpy.asarray(reduced_roots_W, dtype=numpy.complex128)
    if reduced_roots_W.shape != (6,):
        raise ValueError("reduced_roots_W must have shape (6,).")
    if K == 0.0 or sign * K <= 0.0:
        raise ValueError("K and sign are inconsistent.")

    lambda_values = [0.0]
    histories: list[list[complex]] = [[complex(value)] for value in reduced_roots_W]
    current = reduced_roots_W.copy()
    step = float(initial_step)
    rejected = 0
    pseudo_uses = 0

    while lambda_values[-1] < 1.0 - 1.0e-14:
        current_lam = lambda_values[-1]
        target_lam = min(1.0, current_lam + step)
        predictions = numpy.empty(6, dtype=numpy.complex128)
        pseudo_valid = numpy.zeros(6, dtype=bool)

        for branch in range(6):
            previous_lam = lambda_values[-2] if len(lambda_values) >= 2 else None
            previous_W = histories[branch][-2] if len(histories[branch]) >= 2 else None
            seed, valid = _pseudoarclength_homotopy_seed(
                K,
                current_lam,
                histories[branch][-1],
                previous_lam,
                previous_W,
                target_lam,
                parameters,
            )
            predictions[branch] = seed
            pseudo_valid[branch] = valid
        pseudo_uses += int(numpy.count_nonzero(pseudo_valid))

        separation = max(1.0e-8, 1.0e-6 * max(1.0, float(numpy.max(numpy.abs(current)))))
        candidates: list[LocalRoot] = []
        for branch in range(6):
            seed_bank = [predictions[branch], current[branch]]
            if len(lambda_values) >= 2:
                delta_lam = lambda_values[-1] - lambda_values[-2]
                if delta_lam != 0.0:
                    secant = current[branch] + (
                        target_lam - current_lam
                    ) * (current[branch] - histories[branch][-2]) / delta_lam
                    seed_bank.append(secant)
            local_scale = max(0.02, abs(current[branch]), abs(reduced_roots_W[branch]))
            seed_bank.extend(
                predictions[branch] + offset for offset in _broad_seed_offsets(local_scale)
            )
            for seed in seed_bank:
                candidate = solve_homotopy_root_at_lambda(
                    K,
                    target_lam,
                    seed,
                    parameters,
                    residual_tolerance=residual_tolerance,
                )
                _append_unique_candidate(candidates, candidate, separation)

        accepted = False
        if len(candidates) >= 6:
            try:
                assigned, residual, _, _ = _assign_candidates_strict_groups(
                    candidates,
                    predictions,
                    reduced_roots_W,
                    sign,
                )
                accepted = bool(
                    numpy.all(numpy.isfinite(residual))
                    and numpy.max(residual) <= 5.0 * residual_tolerance
                )
            except RuntimeError:
                accepted = False

        if not accepted:
            step *= 0.5
            rejected += 1
            if step < minimum_step:
                raise RuntimeError(
                    f"Homotopy failed at K={K:.8g}, lambda={target_lam:.8g}; "
                    f"only {len(candidates)} distinct roots were found."
                )
            continue

        current = assigned
        lambda_values.append(target_lam)
        for branch in range(6):
            histories[branch].append(complex(current[branch]))
        if numpy.max(residual) < 0.1 * residual_tolerance:
            step = min(maximum_step, 1.35 * step)
        else:
            step = min(maximum_step, 1.12 * step)

    roots, certified, counts = _contour_certify_and_recover(
        K,
        current,
        parameters,
        residual_tolerance=residual_tolerance,
        contour_points=contour_points,
    )
    # Reorder recovered roots against the final homotopy roots.
    candidates = [
        LocalRoot(
            value=complex(value),
            relative_residual=solve_full_root_at_K(
                K, value, parameters, residual_tolerance=residual_tolerance
            ).relative_residual,
            converged=True,
            iterations=0,
        )
        for value in roots
    ]
    roots, residual, _, _ = _assign_candidates_strict_groups(
        candidates, current, reduced_roots_W, sign
    )
    converged = residual <= residual_tolerance

    root_history = numpy.asarray(
        [[histories[branch][index] for branch in range(6)] for index in range(len(lambda_values))],
        dtype=numpy.complex128,
    )
    return HomotopyResult(
        roots=roots,
        residual=residual,
        converged=converged,
        lambda_history=numpy.asarray(lambda_values, dtype=float),
        root_history=root_history,
        pseudoarclength_uses=int(pseudo_uses),
        rejected_steps=int(rejected),
        contour_certified=certified,
        contour_count=counts,
    )


############################################################
# Signed-half-axis continuation
############################################################


def _solve_K_point(
    target_index: int,
    history: list[int],
    K_path: numpy.ndarray,
    reduced_W: numpy.ndarray,
    roots_W: numpy.ndarray,
    parameters: MaxwellianParameters,
    *,
    sign: int,
    slope_order_count: int,
    residual_tolerance: float,
    pseudoarc_stride: int,
    contour_stride: int,
    contour_points: int,
) -> tuple[
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    int,
    bool,
    bool,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
]:
    target_K = float(K_path[target_index])
    prediction = _multiscale_history_predict(
        roots_W, K_path, history, target_index, slope_order_count
    )
    # Keep the reduced ancestry as a lower-weight stabilizer, especially where
    # nearby Landau roots create multiple Newton basins.
    prediction = 0.88 * prediction + 0.12 * reduced_W[target_index]

    pseudo_seed = prediction.copy()
    pseudo_used = False
    if len(history) >= 2 and (
        pseudoarc_stride <= 1 or target_index % pseudoarc_stride == 0
    ):
        previous_index = history[-1]
        earlier_index = history[-2]
        for branch in range(6):
            seed, valid = _pseudoarclength_K_seed(
                float(K_path[previous_index]),
                roots_W[previous_index, branch],
                float(K_path[earlier_index]),
                roots_W[earlier_index, branch],
                target_K,
                parameters,
            )
            if valid:
                pseudo_seed[branch] = seed
                pseudo_used = True

    separation = max(
        1.0e-8,
        1.0e-6 * max(1.0, float(numpy.max(numpy.abs(prediction)))),
    )
    current_index = history[-1]

    # Fast path: one complete-complex correction per predicted sheet.  Recovery
    # seeds are generated only if these six corrections are not distinct.
    primary: list[LocalRoot] = []
    for branch in range(6):
        seed = pseudo_seed[branch] if pseudo_used else prediction[branch]
        candidate = solve_full_root_at_K(
            target_K,
            seed,
            parameters,
            residual_tolerance=residual_tolerance,
        )
        primary.append(candidate)

    candidates: list[LocalRoot] = []
    for candidate in primary:
        _append_unique_candidate(candidates, candidate, separation)

    if len(candidates) < 6:
        # First recovery layer: previous roots and the six reduced seeds.
        for branch in range(6):
            for seed in (
                prediction[branch],
                roots_W[current_index, branch],
                reduced_W[target_index, branch],
            ):
                candidate = solve_full_root_at_K(
                    target_K,
                    seed,
                    parameters,
                    residual_tolerance=residual_tolerance,
                )
                _append_unique_candidate(candidates, candidate, separation)

    if len(candidates) < 6:
        # Second recovery layer: small complex perturbations around only the
        # missing predicted basins.
        for branch in range(6):
            local_scale = max(0.02, abs(prediction[branch]))
            for offset in _broad_seed_offsets(local_scale):
                candidate = solve_full_root_at_K(
                    target_K,
                    prediction[branch] + offset,
                    parameters,
                    residual_tolerance=residual_tolerance,
                )
                _append_unique_candidate(candidates, candidate, separation)
            if len(candidates) >= 6:
                break

    if len(candidates) < 6:
        # Final recovery layer near the three Doppler responses.
        for seed in (
            target_K * drift + imaginary
            for drift in parameters.drift_over_c
            for imaginary in (0.0j, +0.02j, -0.02j, +0.08j, -0.08j)
        ):
            candidate = solve_full_root_at_K(
                target_K,
                seed,
                parameters,
                residual_tolerance=residual_tolerance,
            )
            _append_unique_candidate(candidates, candidate, separation)
            if len(candidates) >= 6:
                break

    assigned, residual, iterations, ties = _assign_candidates_strict_groups(
        candidates,
        prediction,
        reduced_W[target_index],
        sign,
    )

    use_contour = bool(
        contour_stride <= 1
        or target_index % contour_stride == 0
        or numpy.max(residual) > residual_tolerance
        or (
            min(
                abs(assigned[first] - assigned[second])
                for first, second in GROUP_COLUMNS
            ) < 0.015
            and target_index % max(1, contour_stride // 5) == 0
        )
    )
    certified = numpy.zeros(6, dtype=bool)
    counts = numpy.full(6, numpy.nan)
    if use_contour:
        recovered, certified, counts = _contour_certify_and_recover(
            target_K,
            assigned,
            parameters,
            residual_tolerance=residual_tolerance,
            contour_points=contour_points,
        )
        recovery_candidates = []
        for value in recovered:
            item = solve_full_root_at_K(
                target_K,
                value,
                parameters,
                residual_tolerance=residual_tolerance,
            )
            recovery_candidates.append(item)
        if all(item.converged for item in recovery_candidates):
            assigned, residual, iterations, recovery_ties = _assign_candidates_strict_groups(
                recovery_candidates,
                assigned,
                reduced_W[target_index],
                sign,
            )
            ties |= recovery_ties

    return (
        assigned,
        residual,
        iterations,
        len(candidates),
        pseudo_used,
        use_contour,
        certified,
        counts,
        ties,
    )


def solve_full_half_axis(
    k_path: numpy.ndarray,
    reduced_seed: numpy.ndarray,
    parameters: MaxwellianParameters,
    *,
    sign: int,
    anchor_abs_K: float,
    slope_order_count: int,
    residual_tolerance: float,
    pseudoarc_stride: int,
    contour_stride: int,
    contour_points: int,
) -> HalfAxisFullSolution:
    """Solve one open signed half-axis from an independently homotopied anchor."""

    k_path = numpy.asarray(k_path, dtype=float)
    reduced_seed = numpy.asarray(reduced_seed, dtype=numpy.complex128)
    if k_path.ndim != 1 or k_path.size == 0:
        raise ValueError("k_path must be a nonempty one-dimensional array.")
    if reduced_seed.shape != (k_path.size, 6):
        raise ValueError("reduced_seed must have shape (len(k_path),6).")
    if sign not in (-1, +1) or not numpy.all(sign * k_path > 0.0):
        raise ValueError("sign and k_path are inconsistent.")
    if numpy.any(numpy.diff(numpy.abs(k_path)) <= 0.0):
        raise ValueError("k_path must proceed outward from k=0.")

    K_path = k_path * CGS["c"] / parameters.omega_scale
    reduced_W = reduced_seed / parameters.omega_scale
    n_k = k_path.size
    resolved_anchor = min(anchor_abs_K, float(numpy.max(numpy.abs(K_path))))
    anchor_index = int(numpy.argmin(numpy.abs(numpy.abs(K_path) - resolved_anchor)))

    homotopy = solve_anchor_homotopy(
        float(K_path[anchor_index]),
        reduced_W[anchor_index],
        parameters,
        sign=sign,
        residual_tolerance=residual_tolerance,
        contour_points=contour_points,
    )

    roots_W = numpy.full((n_k, 6), numpy.nan + 1j * numpy.nan)
    residual = numpy.full((n_k, 6), numpy.nan)
    converged = numpy.zeros((n_k, 6), dtype=bool)
    iterations = numpy.zeros((n_k, 6), dtype=int)
    candidate_count = numpy.zeros(n_k, dtype=int)
    pseudo_used = numpy.zeros(n_k, dtype=bool)
    contour_used = numpy.zeros(n_k, dtype=bool)
    contour_certified = numpy.zeros((n_k, 6), dtype=bool)
    contour_count = numpy.full((n_k, 6), numpy.nan)
    tie_used = numpy.zeros((n_k, 3), dtype=bool)

    roots_W[anchor_index] = homotopy.roots
    residual[anchor_index] = homotopy.residual
    converged[anchor_index] = homotopy.converged
    iterations[anchor_index] = 0
    candidate_count[anchor_index] = 6
    pseudo_used[anchor_index] = homotopy.pseudoarclength_uses > 0
    contour_used[anchor_index] = True
    contour_certified[anchor_index] = homotopy.contour_certified
    contour_count[anchor_index] = homotopy.contour_count

    history = [anchor_index]
    for index in range(anchor_index + 1, n_k):
        (
            roots_W[index],
            residual[index],
            iterations[index],
            candidate_count[index],
            pseudo_used[index],
            contour_used[index],
            contour_certified[index],
            contour_count[index],
            tie_used[index],
        ) = _solve_K_point(
            index,
            history,
            K_path,
            reduced_W,
            roots_W,
            parameters,
            sign=sign,
            slope_order_count=slope_order_count,
            residual_tolerance=residual_tolerance,
            pseudoarc_stride=pseudoarc_stride,
            contour_stride=contour_stride,
            contour_points=contour_points,
        )
        converged[index] = residual[index] <= residual_tolerance
        history.append(index)

    history = [anchor_index]
    for index in range(anchor_index - 1, -1, -1):
        (
            roots_W[index],
            residual[index],
            iterations[index],
            candidate_count[index],
            pseudo_used[index],
            contour_used[index],
            contour_certified[index],
            contour_count[index],
            tie_used[index],
        ) = _solve_K_point(
            index,
            history,
            K_path,
            reduced_W,
            roots_W,
            parameters,
            sign=sign,
            slope_order_count=slope_order_count,
            residual_tolerance=residual_tolerance,
            pseudoarc_stride=pseudoarc_stride,
            contour_stride=contour_stride,
            contour_points=contour_points,
        )
        converged[index] = residual[index] <= residual_tolerance
        history.append(index)

    return HalfAxisFullSolution(
        k_path=k_path.copy(),
        roots=roots_W * parameters.omega_scale,
        reduced_seed=reduced_seed.copy(),
        relative_residual=residual,
        converged=converged,
        iterations=iterations,
        candidate_count=candidate_count,
        pseudoarclength_used=pseudo_used,
        contour_used=contour_used,
        contour_certified=contour_certified,
        contour_count=contour_count,
        analytic_tie_break_used=tie_used,
        anchor_index=anchor_index,
        homotopy=homotopy,
    )


############################################################
# Full-grid organization and diagnostics
############################################################


def _constant_half_axis_group_relabel(
    negative_aligned: numpy.ndarray,
    positive: numpy.ndarray,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Choose one constant within-group permutation on the negative half-axis."""

    result = negative_aligned.copy()
    swaps = numpy.zeros(3, dtype=bool)
    partner = numpy.asarray(OPPOSITE_K_PARTNER, dtype=int)
    target = -numpy.conjugate(positive[:, partner])

    for group_index, columns in enumerate(GROUP_COLUMNS):
        first, second = columns
        direct = numpy.sum(
            numpy.abs(result[:, first] - target[:, first])
            + numpy.abs(result[:, second] - target[:, second])
        )
        crossed = numpy.sum(
            numpy.abs(result[:, second] - target[:, first])
            + numpy.abs(result[:, first] - target[:, second])
        )
        if crossed + 1.0e-12 < direct:
            result[:, [first, second]] = result[:, [second, first]]
            swaps[group_index] = True
    return result, swaps


def _root_set_mismatch(reference: numpy.ndarray, candidate: numpy.ndarray) -> float:
    maximum = 0.0
    for reference_row, candidate_row in zip(reference, candidate):
        distances = numpy.abs(reference_row[:, None] - candidate_row[None, :])
        mismatch = max(
            float(numpy.max(numpy.min(distances, axis=1))),
            float(numpy.max(numpy.min(distances, axis=0))),
        )
        maximum = max(maximum, mismatch)
    return maximum


def _pair_diagnostics(
    roots: numpy.ndarray,
    scale: float,
    *,
    same_real_tolerance: float,
    opposite_imag_tolerance: float,
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    same_real = numpy.full((roots.shape[0], 2), numpy.nan)
    opposite_imag = numpy.full_like(same_real, numpy.nan)
    magnitude_error = numpy.full_like(same_real, numpy.nan)
    mask = numpy.zeros_like(same_real, dtype=bool)

    for group_index, first in enumerate((2, 4)):
        one = roots[:, first]
        two = roots[:, first + 1]
        same_real[:, group_index] = numpy.abs(one.real - two.real) / scale
        opposite_imag[:, group_index] = numpy.abs(one.imag + two.imag) / scale
        magnitude_error[:, group_index] = numpy.abs(
            numpy.abs(one.imag) - numpy.abs(two.imag)
        ) / scale
        mask[:, group_index] = (
            numpy.isfinite(one.real)
            & numpy.isfinite(two.real)
            & (one.imag * two.imag < 0.0)
            & (same_real[:, group_index] <= same_real_tolerance)
            & (opposite_imag[:, group_index] <= opposite_imag_tolerance)
        )
    return same_real, opposite_imag, magnitude_error, mask


def _contiguous_boundary(k_side: numpy.ndarray, mask_side: numpy.ndarray) -> float:
    order = numpy.argsort(numpy.abs(k_side))
    magnitudes = numpy.abs(k_side[order])
    mask = mask_side[order]
    if magnitudes.size == 0 or not mask[0]:
        return 0.0
    false = numpy.flatnonzero(~mask)
    if false.size == 0:
        return numpy.nan
    index = int(false[0])
    if index == 0:
        return 0.0
    return float(0.5 * (magnitudes[index - 1] + magnitudes[index]))


def _matched_boundary(first: float, second: float) -> float:
    if numpy.isfinite(first) and numpy.isfinite(second):
        return min(first, second)
    return numpy.nan


def _rule5_slope_diagnostics(
    k: numpy.ndarray,
    roots: numpy.ndarray,
    boundary_k: numpy.ndarray,
    slope_order_count: int,
) -> tuple[
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    numpy.ndarray,
    float,
]:
    orders = numpy.arange(1, slope_order_count + 1, dtype=int)
    shape = (4, slope_order_count, 6)
    real_left = numpy.full(shape, numpy.nan)
    real_right = numpy.full(shape, numpy.nan)
    imag_left = numpy.full(shape, numpy.nan)
    imag_right = numpy.full(shape, numpy.nan)
    real_mismatch = numpy.full(shape, numpy.nan)
    imag_mismatch = numpy.full(shape, numpy.nan)
    velocity_floor = 1.0e-12 * CGS["c"]

    for boundary_index, requested in enumerate(boundary_k):
        if not numpy.isfinite(requested) or requested == 0.0:
            continue
        center = int(numpy.argmin(numpy.abs(k - requested)))
        for order_index, order in enumerate(orders):
            left_index = center - int(order)
            right_index = center + int(order)
            if left_index < 0 or right_index >= k.size:
                continue
            left_dk = k[center] - k[left_index]
            right_dk = k[right_index] - k[center]
            if left_dk == 0.0 or right_dk == 0.0:
                continue
            left = (roots[center] - roots[left_index]) / left_dk
            right = (roots[right_index] - roots[center]) / right_dk
            real_left[boundary_index, order_index] = left.real
            real_right[boundary_index, order_index] = right.real
            imag_left[boundary_index, order_index] = left.imag
            imag_right[boundary_index, order_index] = right.imag
            real_scale = numpy.maximum.reduce(
                [
                    numpy.abs(left.real),
                    numpy.abs(right.real),
                    numpy.full(6, velocity_floor),
                ]
            )
            imag_scale = numpy.maximum.reduce(
                [
                    numpy.abs(left.imag),
                    numpy.abs(right.imag),
                    numpy.full(6, velocity_floor),
                ]
            )
            real_mismatch[boundary_index, order_index] = (
                numpy.abs(left.real - right.real) / real_scale
            )
            imag_mismatch[boundary_index, order_index] = (
                numpy.abs(left.imag - right.imag) / imag_scale
            )

    finite = numpy.concatenate(
        [
            real_mismatch[numpy.isfinite(real_mismatch)],
            imag_mismatch[numpy.isfinite(imag_mismatch)],
        ]
    )
    maximum = float(numpy.max(finite)) if finite.size else numpy.nan
    return (
        orders,
        real_left,
        real_right,
        imag_left,
        imag_right,
        real_mismatch,
        imag_mismatch,
        maximum,
    )


def _magnitude_order_departure_count(
    roots: numpy.ndarray, scale: float
) -> tuple[int, int]:
    """Return diagnostic departures from decreasing |Re(omega)| order.

    These counts do not modify roots and are not convergence failures.  If the
    absolute real frequencies of two distinct continuous sheets cross away
    from a coalescence, enforcing the inequality pointwise would require an
    index exchange and would make both Re(omega) and Im(omega) discontinuous.
    """

    tolerance = _CANONICAL_REAL_ORDER_TOL * scale
    g2 = numpy.abs(roots[:, 2].real) + tolerance < numpy.abs(roots[:, 3].real)
    g3 = numpy.abs(roots[:, 4].real) + tolerance < numpy.abs(roots[:, 5].real)
    return int(numpy.count_nonzero(g2)), int(numpy.count_nonzero(g3))


def solve_six_full_maxwellian_branches(
    k: numpy.ndarray,
    reduced_solution: SixComplexRoots,
    wpe_negative: float,
    vd_negative: float,
    vthe_negative: float,
    wpe_positive: float,
    vd_positive: float,
    vthe_positive: float,
    wpi: float,
    vthi: float,
    *,
    anchor_abs_K: float = 5.0,
    slope_order_count: int = 4,
    residual_tolerance: float = 5.0e-11,
    pseudoarc_stride: int = 25,
    contour_stride: int = 250,
    contour_points: int = 96,
    pair_same_real_tolerance: float = 2.0e-3,
    pair_opposite_imag_tolerance: float = 2.0e-3,
) -> SixComplexRoots:
    """Solve six selected full-Maxwellian roots on matched signed half-axes."""

    parameters = build_maxwellian_parameters(
        wpe_negative,
        vd_negative,
        vthe_negative,
        wpe_positive,
        vd_positive,
        vthe_positive,
        wpi,
        vthi,
    )
    k = numpy.asarray(k, dtype=float)
    if k.ndim != 1 or k.size < 3 or numpy.any(numpy.diff(k) <= 0.0):
        raise ValueError("k must be a strictly increasing one-dimensional grid.")
    zero_indices = numpy.flatnonzero(k == 0.0)
    if zero_indices.size != 1:
        raise ValueError("k must contain exactly one zero.")
    zero_index = int(zero_indices[0])
    if reduced_solution.k.shape != k.shape or not numpy.array_equal(reduced_solution.k, k):
        raise ValueError("k must exactly match reduced_solution.k.")
    if slope_order_count < 1:
        raise ValueError("slope_order_count must be at least one.")

    negative_indices = numpy.flatnonzero(k < 0.0)
    positive_indices = numpy.flatnonzero(k > 0.0)
    if negative_indices.size != positive_indices.size:
        raise ValueError("The signed half-axes must have equal lengths.")
    if not numpy.allclose(
        -k[negative_indices][::-1], k[positive_indices], rtol=5.0e-13, atol=0.0
    ):
        raise ValueError("The grid must be matched under k -> -k.")

    positive = solve_full_half_axis(
        k[positive_indices],
        reduced_solution.roots[positive_indices],
        parameters,
        sign=+1,
        anchor_abs_K=anchor_abs_K,
        slope_order_count=slope_order_count,
        residual_tolerance=residual_tolerance,
        pseudoarc_stride=pseudoarc_stride,
        contour_stride=contour_stride,
        contour_points=contour_points,
    )
    negative = solve_full_half_axis(
        k[negative_indices][::-1],
        reduced_solution.roots[negative_indices][::-1],
        parameters,
        sign=-1,
        anchor_abs_K=anchor_abs_K,
        slope_order_count=slope_order_count,
        residual_tolerance=residual_tolerance,
        pseudoarc_stride=pseudoarc_stride,
        contour_stride=contour_stride,
        contour_points=contour_points,
    )

    # Align negative rows by increasing positive |k| and choose only a constant
    # within-group relabeling.  This cannot create a pointwise index exchange.
    negative_aligned, negative_group_swaps = _constant_half_axis_group_relabel(
        negative.roots,
        positive.roots,
    )

    # Apply the same constant negative-side group permutations to diagnostics
    # before the reciprocal signed-k correction below.
    negative_residual = negative.relative_residual.copy()
    negative_converged = negative.converged.copy()
    negative_iterations = negative.iterations.copy()
    negative_certified = negative.contour_certified.copy()
    negative_counts = negative.contour_count.copy()
    for group_index, swapped in enumerate(negative_group_swaps):
        if swapped:
            columns = list(GROUP_COLUMNS[group_index])
            negative_residual[:, columns] = negative_residual[:, columns[::-1]]
            negative_converged[:, columns] = negative_converged[:, columns[::-1]]
            negative_iterations[:, columns] = negative_iterations[:, columns[::-1]]
            negative_certified[:, columns] = negative_certified[:, columns[::-1]]
            negative_counts[:, columns] = negative_counts[:, columns[::-1]]

    # The independent half-axis continuations can settle on neighboring damped
    # Landau zeros even though both roots have small residuals.  Reconcile only
    # such signed-k disagreements by using the exact transformed counterpart as
    # an additional seed and correcting it directly at the negative wavenumber.
    # No root is copied: every accepted value is re-solved against the negative-k
    # Maxwellian equation.  This preserves the exact signed-k branch identity
    # while retaining the independent half-axis calculations as the first pass.
    partner = numpy.asarray(OPPOSITE_K_PARTNER, dtype=int)
    reciprocal_tolerance = max(5.0e-10, 20.0 * residual_tolerance)
    for row_index, k_value in enumerate(negative.k_path):
        K_negative = k_value * CGS["c"] / parameters.omega_scale
        for branch in range(6):
            target_W = -numpy.conjugate(
                positive.roots[row_index, partner[branch]]
            ) / parameters.omega_scale
            disagreement = abs(
                negative_aligned[row_index, branch]
                / parameters.omega_scale
                - target_W
            )
            if disagreement <= reciprocal_tolerance:
                continue
            corrected = solve_full_root_at_K(
                K_negative,
                target_W,
                parameters,
                residual_tolerance=residual_tolerance,
            )
            if corrected.converged:
                negative_aligned[row_index, branch] = (
                    corrected.value * parameters.omega_scale
                )
                negative_residual[row_index, branch] = corrected.relative_residual
                negative_converged[row_index, branch] = True
                negative_iterations[row_index, branch] = corrected.iterations

    roots = numpy.full((k.size, 6), numpy.nan + 1j * numpy.nan)
    reduced_seed = reduced_solution.roots.copy()
    residual = numpy.full((k.size, 6), numpy.nan)
    converged = numpy.zeros((k.size, 6), dtype=bool)
    iterations = numpy.zeros((k.size, 6), dtype=int)
    solved_directly = numpy.zeros((k.size, 6), dtype=bool)
    candidate_count = numpy.zeros(k.size, dtype=int)
    pseudo_used = numpy.zeros(k.size, dtype=bool)
    contour_used = numpy.zeros(k.size, dtype=bool)
    contour_certified = numpy.zeros((k.size, 6), dtype=bool)
    contour_count = numpy.full((k.size, 6), numpy.nan)
    tie_used = numpy.zeros((k.size, 3), dtype=bool)

    roots[positive_indices] = positive.roots
    roots[negative_indices] = negative_aligned[::-1]
    total_frequency = float(
        numpy.sqrt(wpe_negative**2 + wpe_positive**2 + wpi**2)
    )
    roots[zero_index] = numpy.asarray(
        [+total_frequency, -total_frequency, 0.0j, 0.0j, 0.0j, 0.0j]
    )

    residual[positive_indices] = positive.relative_residual
    residual[negative_indices] = negative_residual[::-1]
    converged[positive_indices] = positive.converged
    converged[negative_indices] = negative_converged[::-1]
    iterations[positive_indices] = positive.iterations
    iterations[negative_indices] = negative_iterations[::-1]
    solved_directly[positive_indices] = True
    solved_directly[negative_indices] = True
    candidate_count[positive_indices] = positive.candidate_count
    candidate_count[negative_indices] = negative.candidate_count[::-1]
    pseudo_used[positive_indices] = positive.pseudoarclength_used
    pseudo_used[negative_indices] = negative.pseudoarclength_used[::-1]
    contour_used[positive_indices] = positive.contour_used
    contour_used[negative_indices] = negative.contour_used[::-1]
    contour_certified[positive_indices] = positive.contour_certified
    contour_certified[negative_indices] = negative_certified[::-1]
    contour_count[positive_indices] = positive.contour_count
    contour_count[negative_indices] = negative_counts[::-1]
    tie_used[positive_indices] = positive.analytic_tie_break_used
    tie_used[negative_indices] = negative.analytic_tie_break_used[::-1]

    same_real, opposite_imag, magnitude_error, pair_mask = _pair_diagnostics(
        roots,
        parameters.omega_scale,
        same_real_tolerance=pair_same_real_tolerance,
        opposite_imag_tolerance=pair_opposite_imag_tolerance,
    )

    full_k_b1_negative = _contiguous_boundary(
        k[negative_indices], pair_mask[negative_indices, 0]
    )
    full_k_b1_positive = _contiguous_boundary(
        k[positive_indices], pair_mask[positive_indices, 0]
    )
    full_k_b2_negative = _contiguous_boundary(
        k[negative_indices], pair_mask[negative_indices, 1]
    )
    full_k_b2_positive = _contiguous_boundary(
        k[positive_indices], pair_mask[positive_indices, 1]
    )
    full_k_b1_common = _matched_boundary(full_k_b1_negative, full_k_b1_positive)
    full_k_b2_common = _matched_boundary(full_k_b2_negative, full_k_b2_positive)

    k_b1_negative = reduced_solution.k_b1_negative
    k_b1_positive = reduced_solution.k_b1_positive
    k_b1_common = reduced_solution.k_b1_common
    k_b2_negative = reduced_solution.k_b2_negative
    k_b2_positive = reduced_solution.k_b2_positive
    k_b2_common = reduced_solution.k_b2_common

    boundary_k = numpy.asarray(
        [-k_b1_negative, +k_b1_positive, -k_b2_negative, +k_b2_positive],
        dtype=float,
    )
    (
        slope_orders,
        slope_real_left,
        slope_real_right,
        slope_imag_left,
        slope_imag_right,
        slope_real_mismatch,
        slope_imag_mismatch,
        maximum_slope_mismatch,
    ) = _rule5_slope_diagnostics(k, roots, boundary_k, slope_order_count)

    aligned_negative = roots[negative_indices][::-1]
    target = -numpy.conjugate(
        roots[positive_indices][:, numpy.asarray(OPPOSITE_K_PARTNER, dtype=int)]
    )
    opposite_error = numpy.max(numpy.abs(aligned_negative - target), axis=0)
    root_set_mismatch = _root_set_mismatch(target, aligned_negative)

    nonzero = k != 0.0
    finite_bounds = [
        value
        for value in (k_b1_common, k_b2_common)
        if numpy.isfinite(value) and value > 0.0
    ]
    common_limit = min(finite_bounds) if finite_bounds else numpy.inf
    hierarchy_mask = (
        nonzero
        & pair_mask[:, 0]
        & pair_mask[:, 1]
        & (numpy.abs(k) < common_limit)
    )
    g2 = 0.5 * (numpy.abs(roots[:, 2].imag) + numpy.abs(roots[:, 3].imag))
    g3 = 0.5 * (numpy.abs(roots[:, 4].imag) + numpy.abs(roots[:, 5].imag))
    hierarchy_margin = (g2 - g3) / parameters.omega_scale
    if numpy.any(hierarchy_mask):
        minimum_margin = float(numpy.min(hierarchy_margin[hierarchy_mask]))
        violations = int(numpy.count_nonzero(hierarchy_margin[hierarchy_mask] <= 0.0))
    else:
        minimum_margin = numpy.nan
        violations = 0

    # The magnitude-order counts are retained as diagnostics only.  Complete-
    # root continuity has priority over a pointwise |Re(omega)| ranking.
    _magnitude_order_departure_count(roots, parameters.omega_scale)

    k_bs = numpy.asarray(
        [
            [full_k_b1_negative, full_k_b1_positive, full_k_b1_common],
            [full_k_b2_negative, full_k_b2_positive, full_k_b2_common],
        ],
        dtype=float,
    )
    reference_k_bs = numpy.asarray(
        [
            [k_b1_negative, k_b1_positive, k_b1_common],
            [k_b2_negative, k_b2_positive, k_b2_common],
        ],
        dtype=float,
    )
    solution_parameters = {
        "model": "cold_ion_two_electron_maxwellian",
        "wpe_negative": float(wpe_negative),
        "vd_negative": float(vd_negative),
        "vthe_negative": float(vthe_negative),
        "wpe_positive": float(wpe_positive),
        "vd_positive": float(vd_positive),
        "vthe_positive": float(vthe_positive),
        "wpi": float(wpi),
        "vthi": float(vthi),
        "anchor_abs_K": float(anchor_abs_K),
        "slope_order_count": int(slope_order_count),
        "residual_tolerance": float(residual_tolerance),
        "pseudoarc_stride": int(pseudoarc_stride),
        "contour_stride": int(contour_stride),
        "contour_points": int(contour_points),
        "pair_same_real_tolerance": float(pair_same_real_tolerance),
        "pair_opposite_imag_tolerance": float(pair_opposite_imag_tolerance),
    }

    return SixComplexRoots(
        k=k.copy(),
        roots=roots,
        k_bs=k_bs,
        parameters=solution_parameters,
        reference_k_bs=reference_k_bs,
        reduced_seed=reduced_seed,
        relative_residual=residual,
        converged=converged,
        iterations=iterations,
        solved_directly=solved_directly,
        candidate_count=candidate_count,
        pseudoarclength_used=pseudo_used,
        contour_used=contour_used,
        contour_certified=contour_certified,
        contour_count=contour_count,
        analytic_tie_break_used=tie_used,
        positive_anchor_index=int(positive.anchor_index),
        negative_anchor_index=int(negative.anchor_index),
        positive_homotopy=positive.homotopy,
        negative_homotopy=negative.homotopy,
        same_real_error=same_real,
        opposite_imag_error=opposite_imag,
        growth_magnitude_error=magnitude_error,
        pair_condition_mask=pair_mask,
        k_b1_negative=k_b1_negative,
        k_b1_positive=k_b1_positive,
        k_b1_common=k_b1_common,
        k_b2_negative=k_b2_negative,
        k_b2_positive=k_b2_positive,
        k_b2_common=k_b2_common,
        full_k_b1_negative=full_k_b1_negative,
        full_k_b1_positive=full_k_b1_positive,
        full_k_b1_common=full_k_b1_common,
        full_k_b2_negative=full_k_b2_negative,
        full_k_b2_positive=full_k_b2_positive,
        full_k_b2_common=full_k_b2_common,
        growth_hierarchy_mask=hierarchy_mask,
        growth_hierarchy_margin=hierarchy_margin,
        minimum_growth_hierarchy_margin=minimum_margin,
        growth_hierarchy_violation_count=violations,
        opposite_k_conversion_error=opposite_error,
        opposite_k_root_set_mismatch=float(root_set_mismatch),
        slope_boundary_k=boundary_k,
        slope_orders=slope_orders,
        slope_real_left=slope_real_left,
        slope_real_right=slope_real_right,
        slope_imag_left=slope_imag_left,
        slope_imag_right=slope_imag_right,
        slope_real_relative_mismatch=slope_real_mismatch,
        slope_imag_relative_mismatch=slope_imag_mismatch,
        maximum_rule5_relative_mismatch=maximum_slope_mismatch,
    )







############################################################
############################################################

def _hdf5_parameter_values(solution: SixComplexRoots) -> dict[str, float]:
    parameters = solution.parameters
    wpe1 = float(parameters["wpe_negative"])
    wpe2 = float(parameters["wpe_positive"])
    wpe = float(numpy.hypot(wpe1, wpe2))

    return {
        "wpe": wpe,
        "wpe1": wpe1,
        "wpe2": wpe2,
        "vd1": float(parameters["vd_negative"]),
        "vthe1": float(parameters.get("vthe_negative", 0.0)),
        "vd2": float(parameters["vd_positive"]),
        "vthe2": float(parameters.get("vthe_positive", 0.0)),
        "wpi": float(parameters["wpi"]),
        "vthi": float(parameters.get("vthi", 0.0)),
    }


def save_solution_hdf5(
    path: str | Path,
    solution: SixComplexRoots,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    reference_kb = solution.reference_k_bs
    if reference_kb is None:
        reference_kb = solution.k_bs

    with h5py.File(path, "w") as stream:
        parameter_group = stream.create_group("parameters", track_order=True)
        for name, value in _hdf5_parameter_values(solution).items():
            parameter_group.create_dataset(name, data=value)

        kb = numpy.array([solution.k_b1_common,solution.k_b2_common],dtype=float)

        root_group = stream.create_group("roots", track_order=True)
        root_group.create_dataset("k", data=solution.k)
        root_group.create_dataset("kb", data=kb)
        root_group.create_dataset("imag", data=solution.roots_imag, compression="gzip", shuffle=True)
        root_group.create_dataset("real", data=solution.roots_real, compression="gzip", shuffle=True)






############################################################
# Bohm-Gross Langmuir mode
#
# wave equation : $\omega_L(k)=kv_{b1}+\sqrt{\omega_{pe1}^2+3k^2v_{\mathrm{the1}}^2}$
#
############################################################
def dispersion_relation_langmuir(k,wpe1,vd1,vthe1):
    return k*vd1+numpy.sqrt(wpe1**2+3.0*numpy.power(k*vthe1,2.0))



def growth_rate_langmuir(k,wpe1,vd1,vthe1):
    wL   = dispersion_relation_langmuir(k,wpe1,vd1,vthe1)

    NL   = wpe1**2*(2.0*numpy.power(wL-k*vd1,-3.0)+12.0*numpy.power(k*vthe1,2.0)*numpy.power(wL-k*vd1,-5.0))

    xiL1 = (wL-k*vd1)/vthe1/numpy.abs(k)/numpy.sqrt(2)

    IL   = numpy.sqrt(numpy.pi)*(numpy.power(wpe1/k/vthe1,2.0)*xiL1*numpy.exp(-xiL1**2))

    return -1.0*IL/NL





############################################################
# Ion-acoustic mode
#
# wave equation : $\mathcal{D}_{R,\mathrm{IA}}(k,\omega_r) = 1 - \frac{\omega_{pi}^{2}}{\omega_r^{2}} \left(1+\frac{3k^{2}v_{\mathrm{thi}}^{2}}{\omega_r^{2}}\right) + \frac{1}{k^{2}\lambda_D^{2}}$
#
############################################################
def dispersion_relation_ion_acoustic(k,wpe1,vthe1,wpi):

    lambdaD1 = vthe1/wpe1
    
    cs       = wpi*lambdaD1

    return cs,1.0*cs*numpy.abs(k)/numpy.sqrt(1.0+numpy.power(k*lambdaD1,2.0))


def growth_rate_ion_acoustic(k,wpe1,vd1,vthe1,wpi,vthi):
    
    _,wIA = dispersion_relation_ion_acoustic(k,wpe1,vthe1,wpi)


    NIA   = 2.0*numpy.power(wpi,2.0)*numpy.power(wIA,-3.0)+12.0*numpy.power(wpi*k*vthi,2.0)*numpy.power(wIA,-5.0)

    xii   = wIA/vthi/numpy.abs(k)/numpy.sqrt(2)
    xie1  = (wIA-k*vd1)/vthe1/numpy.abs(k)/numpy.sqrt(2)

    IIA   = numpy.sqrt(numpy.pi)*(numpy.power(wpi/k/vthi,2.0)*xii*numpy.exp(-xii**2)+numpy.power(wpe1/k/vthe1,2.0)*xie1*numpy.exp(-xie1**2))

    return -1.0*IIA/NIA




############################################################
# MHD waves
############################################################

def dispersion_relation_MHD_transverse(k,wpe,c=CGS["c"]):
    return numpy.sqrt(wpe**2+numpy.power(c*k,2.0))
