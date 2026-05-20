import os,sys
sys.path.append(os.getenv("PATH_YAOXPY"))



############################################################
def fun_wavenumber_sym(n,dk):
    return numpy.arange(-n,n+1)*dk



############################################################
import numpy
from scipy.optimize import linear_sum_assignment


# ============================================================
# shared small utilities
# ============================================================
_EPS = 1.0e-30


def _finite_diff_complex(y0, y1, x0, x1):
    dx = x1 - x0
    if abs(dx) < _EPS:
        return 0.0 + 0.0j
    return (y1 - y0) / dx


def _finite_diff_real(y0, y1, x0, x1):
    dx = x1 - x0
    if abs(dx) < _EPS:
        return 0.0
    return (y1 - y0) / dx


def _linear_predict(y2, y1, x2, x1, x0):
    dx = x1 - x2
    if abs(dx) < _EPS:
        return y1
    slope = (y1 - y2) / dx
    return y1 + slope * (x0 - x1)


def _sort_initial_roots(w0, gfun, a=1.0, b=1.0):
    g0 = gfun(w0, a=a, b=b)
    key = numpy.lexsort((numpy.imag(w0), numpy.real(w0), g0))
    return w0[key]


def _choose_start_side(w_all):
    left_val = numpy.max(numpy.abs(numpy.real(w_all[0])))
    right_val = numpy.max(numpy.abs(numpy.real(w_all[-1])))
    return "left" if left_val > right_val else "right"


def _keyfun(x, tol):
    return int(numpy.round(x / tol))


# ============================================================
# generic tracking machinery
# ============================================================
def _branch_transition_cost(
    x2,
    x1,
    x0,
    w2,
    w1,
    w0,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
):
    wp = _linear_predict(w2, w1, x2, x1, x0)

    gp2 = gfun(w2, a=a, b=b)
    gp1 = gfun(w1, a=a, b=b)
    gp0 = gfun(w0, a=a, b=b)
    gp = _linear_predict(gp2, gp1, x2, x1, x0)

    c_dist = abs(w0 - wp) ** 2

    s_prev = _finite_diff_complex(w2, w1, x2, x1)
    s_new = _finite_diff_complex(w1, w0, x1, x0)
    c_slope = abs(s_new - s_prev) ** 2

    c_g = (gp0 - gp) ** 2

    gs_prev = _finite_diff_real(gp2, gp1, x2, x1)
    gs_new = _finite_diff_real(gp1, gp0, x1, x0)
    c_gslope = (gs_new - gs_prev) ** 2

    c_curv = abs(w0 - 2.0 * w1 + w2) ** 2 + (gp0 - 2.0 * gp1 + gp2) ** 2

    return (
        lam_dist * c_dist
        + lam_slope * c_slope
        + lam_g * c_g
        + lam_gslope * c_gslope
        + lam_curv * c_curv
    )


def _build_cost_matrix(
    x2,
    x1,
    x0,
    w2,
    w1,
    wcand,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
):
    nm = len(w1)
    C = numpy.empty((nm, nm), dtype=float)

    for m in range(nm):
        wm = wcand[m]
        for n in range(nm):
            C[m, n] = _branch_transition_cost(
                x2,
                x1,
                x0,
                w2[n],
                w1[n],
                wm,
                gfun,
                a=a,
                b=b,
                lam_dist=lam_dist,
                lam_slope=lam_slope,
                lam_g=lam_g,
                lam_gslope=lam_gslope,
                lam_curv=lam_curv,
            )

    return C


def _track_one_side(
    k_half,
    w_half,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
):
    nk, nm = w_half.shape
    out = numpy.empty_like(w_half, dtype=complex)

    out[0] = _sort_initial_roots(w_half[0], gfun, a=a, b=b)

    if nk == 1:
        return out

    C1 = numpy.empty((nm, nm), dtype=float)
    g0 = gfun(out[0], a=a, b=b)
    g1 = gfun(w_half[1], a=a, b=b)

    for m in range(nm):
        for n in range(nm):
            C1[m, n] = abs(w_half[1, m] - out[0, n]) ** 2 + 0.8 * (g1[m] - g0[n]) ** 2

    row, col = linear_sum_assignment(C1)
    out[1, col] = w_half[1, row]

    for i in range(2, nk):
        C = _build_cost_matrix(
            k_half[i - 2],
            k_half[i - 1],
            k_half[i],
            out[i - 2],
            out[i - 1],
            w_half[i],
            gfun,
            a=a,
            b=b,
            lam_dist=lam_dist,
            lam_slope=lam_slope,
            lam_g=lam_g,
            lam_gslope=lam_gslope,
            lam_curv=lam_curv,
        )
        row, col = linear_sum_assignment(C)
        out[i, col] = w_half[i, row]

    return out


def _pair_step_cost(
    x2,
    x1,
    x0,
    up2,
    up1,
    up0,
    uq2,
    uq1,
    uq0,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
):
    cp = _branch_transition_cost(
        x2,
        x1,
        x0,
        up2,
        up1,
        up0,
        gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
    )

    cq = _branch_transition_cost(
        x2,
        x1,
        x0,
        uq2,
        uq1,
        uq0,
        gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
    )

    sep_prev = abs(up1 - uq1)
    sep_now = abs(up0 - uq0)
    csep = (sep_now - sep_prev) ** 2

    return cp + cq + lam_sep * csep


def _repair_self_pairs_forward(
    k_half,
    tracked_half,
    self_pairs,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
):
    out = tracked_half.copy()
    nk = len(k_half)

    if nk < 3:
        return out, False

    changed = False

    for i in range(2, nk):
        x2, x1, x0 = k_half[i - 2], k_half[i - 1], k_half[i]

        for p, q in self_pairs:
            keep = _pair_step_cost(
                x2,
                x1,
                x0,
                out[i - 2, p],
                out[i - 1, p],
                out[i, p],
                out[i - 2, q],
                out[i - 1, q],
                out[i, q],
                gfun,
                a=a,
                b=b,
                lam_dist=lam_dist,
                lam_slope=lam_slope,
                lam_g=lam_g,
                lam_gslope=lam_gslope,
                lam_curv=lam_curv,
                lam_sep=lam_sep,
            )

            swap = _pair_step_cost(
                x2,
                x1,
                x0,
                out[i - 2, p],
                out[i - 1, p],
                out[i, q],
                out[i - 2, q],
                out[i - 1, q],
                out[i, p],
                gfun,
                a=a,
                b=b,
                lam_dist=lam_dist,
                lam_slope=lam_slope,
                lam_g=lam_g,
                lam_gslope=lam_gslope,
                lam_curv=lam_curv,
                lam_sep=lam_sep,
            )

            if swap + 1.0e-14 < keep:
                out[i, p], out[i, q] = out[i, q], out[i, p]
                changed = True

    return out, changed


def _repair_self_odd_pairs(
    k_half,
    tracked_half,
    self_pairs,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    passes=8,
):
    out = tracked_half.copy()

    for _ in range(passes):
        changed_any = False

        out, changed = _repair_self_pairs_forward(
            k_half,
            out,
            self_pairs,
            gfun,
            a=a,
            b=b,
            lam_dist=lam_dist,
            lam_slope=lam_slope,
            lam_g=lam_g,
            lam_gslope=lam_gslope,
            lam_curv=lam_curv,
            lam_sep=lam_sep,
        )
        changed_any = changed_any or changed

        krev = k_half[::-1].copy()
        wrev = out[::-1].copy()

        wrev, changed = _repair_self_pairs_forward(
            krev,
            wrev,
            self_pairs,
            gfun,
            a=a,
            b=b,
            lam_dist=lam_dist,
            lam_slope=lam_slope,
            lam_g=lam_g,
            lam_gslope=lam_gslope,
            lam_curv=lam_curv,
            lam_sep=lam_sep,
        )
        out = wrev[::-1].copy()
        changed_any = changed_any or changed

        if not changed_any:
            break

    return out


def _reconstruct_full_domain_exact(
    k,
    idx_half,
    tracked_half,
    nm,
    cross_pair,
    self_pairs,
    tol=1.0e-12,
):
    k = numpy.asarray(k, dtype=float)
    tracked = numpy.zeros((len(k), nm), dtype=complex)
    tracked[idx_half] = tracked_half

    idx_map = {_keyfun(ki, tol): i for i, ki in enumerate(k)}
    half_set = set(idx_half.tolist())

    c0, c1 = cross_pair

    for i in idx_half:
        km = -k[i]
        kk = _keyfun(km, tol)
        if kk not in idx_map:
            continue

        j = idx_map[kk]
        if j in half_set:
            continue

        tracked[j, c0] = -tracked[i, c1]
        tracked[j, c1] = -tracked[i, c0]

        for p, q in self_pairs:
            tracked[j, p] = -tracked[i, p]
            tracked[j, q] = -tracked[i, q]

    return tracked


def _track_branches_symmetric(
    k,
    w_all,
    nm,
    cross_pair,
    self_pairs,
    gfun,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8,
):
    k = numpy.asarray(k, dtype=float)
    w_all = numpy.asarray(w_all, dtype=complex)

    if w_all.shape != (len(k), nm):
        raise ValueError(f"w_all must have shape (len(k), {nm})")

    side = _choose_start_side(w_all)

    if side == "right":
        idx_half = numpy.where(k >= -tol)[0]
        idx_half = idx_half[numpy.argsort(k[idx_half])[::-1]]
    else:
        idx_half = numpy.where(k <= tol)[0]
        idx_half = idx_half[numpy.argsort(k[idx_half])]

    k_half = k[idx_half]
    w_half = w_all[idx_half]

    tracked_half = _track_one_side(
        k_half,
        w_half,
        gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
    )

    tracked_half = _repair_self_odd_pairs(
        k_half,
        tracked_half,
        self_pairs,
        gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
        lam_sep=lam_sep,
        passes=repair_passes,
    )

    tracked_full = _reconstruct_full_domain_exact(
        k,
        idx_half,
        tracked_half,
        nm,
        cross_pair,
        self_pairs,
        tol=tol,
    )

    return tracked_full.T


def _check_symmetry(k, omega, cross_pair, self_pairs, gfun, a=1.0, b=1.0, tol=1.0e-12):
    k = numpy.asarray(k, dtype=float)
    omega = numpy.asarray(omega, dtype=complex)

    idx_map = {_keyfun(ki, tol): i for i, ki in enumerate(k)}
    report = {}

    c0, c1 = cross_pair

    err_cross_complex = []
    err_cross_g = []
    for i, ki in enumerate(k):
        kk = _keyfun(-ki, tol)
        if kk in idx_map:
            j = idx_map[kk]
            err_cross_complex.append(abs(omega[c0, i] + omega[c1, j]))
            err_cross_g.append(abs(gfun(omega[c0, i], a=a, b=b) + gfun(omega[c1, j], a=a, b=b)))

    report[f"pair_{c0}_{c1}_complex_cross_error_max"] = numpy.max(err_cross_complex) if err_cross_complex else numpy.nan
    report[f"pair_{c0}_{c1}_g_cross_error_max"] = numpy.max(err_cross_g) if err_cross_g else numpy.nan

    for p, q in self_pairs:
        errp_complex = []
        errq_complex = []
        errp_g = []
        errq_g = []

        for i, ki in enumerate(k):
            kk = _keyfun(-ki, tol)
            if kk in idx_map:
                j = idx_map[kk]
                errp_complex.append(abs(omega[p, i] + omega[p, j]))
                errq_complex.append(abs(omega[q, i] + omega[q, j]))
                errp_g.append(abs(gfun(omega[p, i], a=a, b=b) + gfun(omega[p, j], a=a, b=b)))
                errq_g.append(abs(gfun(omega[q, i], a=a, b=b) + gfun(omega[q, j], a=a, b=b)))

        report[f"pair_{p}_{q}_branch_{p}_complex_self_error_max"] = numpy.max(errp_complex) if errp_complex else numpy.nan
        report[f"pair_{p}_{q}_branch_{q}_complex_self_error_max"] = numpy.max(errq_complex) if errq_complex else numpy.nan
        report[f"pair_{p}_{q}_branch_{p}_g_self_error_max"] = numpy.max(errp_g) if errp_g else numpy.nan
        report[f"pair_{p}_{q}_branch_{q}_g_self_error_max"] = numpy.max(errq_g) if errq_g else numpy.nan

    return report


def _check_local_slope_jumps(k, omega):
    k = numpy.asarray(k, dtype=float)
    omega = numpy.asarray(omega, dtype=complex)

    out = {}
    nk = len(k)

    if nk < 3:
        for n in range(omega.shape[0]):
            out[f"branch_{n}_max_second_diff"] = numpy.nan
        return out

    for n in range(omega.shape[0]):
        vals = []
        for i in range(2, nk):
            vals.append(abs(omega[n, i] - 2.0 * omega[n, i - 1] + omega[n, i - 2]))
        out[f"branch_{n}_max_second_diff"] = max(vals)

    return out


# ============================================================
# 12th-order coefficients
# note: left mostly algebraically unchanged for safety
# ============================================================
def wave_equation_two_electrons_twelveth_coefficients(wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k):
    wpi = wpe / numpy.sqrt(mu)
    vthi0 = vthe0 / numpy.sqrt(mu)
    vthi1 = vthe1 / numpy.sqrt(mu)

    a12 = -1.0

    a11 = 4.0 * (vd0 + vd1) * k

    a10 = -(6.0 * vd1 * vd1 + 16.0 * vd0 * vd1 + 6.0 * vd0 * vd0) * k * k \
          + (alpha0 + alpha1) * (wpi * wpi + wpe * wpe)

    a09 = 4.0 * (
        vd1 * vd1 * vd1 + 6.0 * vd0 * vd1 * vd1 + 6.0 * vd0 * vd0 * vd1 + vd0 * vd0 * vd0
    ) * k * k * k \
          - 4.0 * wpi * wpi * (alpha0 + alpha1) * (vd0 + vd1) * k \
          - 2.0 * wpe * wpe * (alpha0 * (2.0 * vd1 + vd0) + alpha1 * (vd1 + 2.0 * vd0)) * k

    a08 = -(
        vd1 ** 4 + 16.0 * vd0 * vd1 ** 3 + 36.0 * vd0 * vd0 * vd1 * vd1
        + 16.0 * vd0 ** 3 * vd1 + vd0 ** 4
    ) * k ** 4 \
          + wpi * wpi * (
              (alpha0 + alpha1) * (6.0 * vd1 * vd1 + 16.0 * vd0 * vd1 + 6.0 * vd0 * vd0)
              + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1)
          ) * k * k \
          + wpe * wpe * (
              alpha0 * (6.0 * vd1 * vd1 + 8.0 * vd0 * vd1 + vd0 * vd0) + 3.0 * alpha0 * vthe0 * vthe0
              + alpha1 * (vd1 * vd1 + 8.0 * vd0 * vd1 + 6.0 * vd0 * vd0) + 3.0 * alpha1 * vthe1 * vthe1
          ) * k * k

    a07 = 4.0 * (
        vd0 * vd1 ** 4 + 6.0 * vd0 * vd0 * vd1 ** 3
        + 6.0 * vd0 ** 3 * vd1 * vd1 + vd0 ** 4 * vd1
    ) * k ** 5 \
          - 4.0 * wpi * wpi * (
              (alpha0 + alpha1) * (
                  vd1 ** 3 + 6.0 * vd0 * vd1 * vd1 + 6.0 * vd0 * vd0 * vd1 + vd0 ** 3
              )
              + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1) * (vd0 + vd1)
          ) * k ** 3 \
          - 4.0 * wpe * wpe * (
              alpha0 * (vd1 ** 3 + 3.0 * vd0 * vd1 * vd1 + vd0 * vd0 * vd1)
              + 3.0 * alpha0 * vthe0 * vthe0 * vd1
              + alpha1 * (vd0 * vd1 * vd1 + 3.0 * vd0 * vd0 * vd1 + vd0 ** 3)
              + 3.0 * alpha1 * vthe1 * vthe1 * vd0
          ) * k ** 3

    a06 = -(
        6.0 * vd0 * vd0 * vd1 ** 4 + 16.0 * (vd0 * vd1) ** 3
        + 6.0 * vd0 ** 4 * vd1 * vd1
    ) * k ** 6 \
          + wpi * wpi * (
              (alpha0 + alpha1) * (
                  vd1 ** 4 + 16.0 * vd0 * vd1 ** 3 + 36.0 * (vd0 * vd1) ** 2
                  + 16.0 * vd0 ** 3 * vd1 + vd0 ** 4
              )
              + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1)
              * (6.0 * vd1 * vd1 + 16.0 * vd0 * vd1 + 6.0 * vd0 * vd0)
          ) * k ** 4 \
          + wpe * wpe * (
              alpha0 * (vd1 ** 4 + 8.0 * vd0 * vd1 ** 3 + 6.0 * (vd0 * vd1) ** 2)
              + 18.0 * alpha0 * (vthe0 * vd1) ** 2
              + alpha1 * (6.0 * (vd0 * vd1) ** 2 + 8.0 * vd0 ** 3 * vd1 + vd0 ** 4)
              + 18.0 * alpha1 * (vthe1 * vd0) ** 2
          ) * k ** 4

    a05 = 4.0 * (
        vd0 ** 3 * vd1 ** 4 + vd0 ** 4 * vd1 ** 3
    ) * k ** 7 \
          - 4.0 * wpi * wpi * (
              (alpha0 + alpha1) * (
                  vd0 * vd1 ** 4 + 6.0 * vd0 * vd0 * vd1 ** 3
                  + 6.0 * vd0 ** 3 * vd1 * vd1 + vd0 ** 4 * vd1
              )
              + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1)
              * (vd1 ** 3 + 6.0 * vd0 * vd1 * vd1 + 6.0 * vd0 * vd0 * vd1 + vd0 ** 3)
          ) * k ** 5 \
          - 2.0 * wpe * wpe * (
              alpha0 * (vd0 * vd1 ** 4 + 2.0 * vd0 * vd0 * vd1 ** 3)
              + 6.0 * alpha0 * vthe0 * vthe0 * vd1 ** 3
              + alpha1 * (2.0 * vd0 ** 3 * vd1 * vd1 + vd0 ** 4 * vd1)
              + 6.0 * alpha1 * vthe1 * vthe1 * vd0 ** 3
          ) * k ** 5

    a04 = -(vd0 * vd1) ** 4 * k ** 8 \
          + wpi * wpi * (
              (alpha0 + alpha1) * (
                  6.0 * (vd0 * vd1 * vd1) ** 2 + 16.0 * (vd0 * vd1) ** 3
                  + 6.0 * (vd0 * vd0 * vd1) ** 2
              )
              + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1)
              * (
                  vd1 ** 4 + 16.0 * vd0 * vd1 ** 3 + 36.0 * (vd0 * vd1) ** 2
                  + 16.0 * vd0 ** 3 * vd1 + vd0 ** 4
              )
          ) * k ** 6 \
          + wpe * wpe * (
              alpha0 * (vd0 * vd1 * vd1) ** 2 + 3.0 * alpha0 * (vthe0 * vd1 * vd1) ** 2
              + alpha1 * (vd0 * vd0 * vd1) ** 2 + 3.0 * alpha1 * (vthe1 * vd0 * vd0) ** 2
          ) * k ** 6

    a03 = -4.0 * wpi * wpi * (
        (alpha0 + alpha1) * (
            vd0 ** 3 * vd1 ** 4 + vd0 ** 4 * vd1 ** 3
        )
        + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1) * (
            vd0 * vd1 ** 4 + 6.0 * vd0 * vd0 * vd1 ** 3
            + 6.0 * vd0 ** 3 * vd1 * vd1 + vd0 ** 4 * vd1
        )
    ) * k ** 7

    a02 = wpi * wpi * (
        (alpha0 + alpha1) * (vd0 * vd1) ** 4
        + 3.0 * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1)
        * (
            6.0 * (vd0 * vd1 * vd1) ** 2 + 16.0 * (vd0 * vd1) ** 3
            + 6.0 * (vd0 * vd0 * vd1) ** 2
        )
    ) * k ** 8

    a01 = -12.0 * wpi * wpi * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1) * (
        vd0 ** 3 * vd1 ** 4 + vd0 ** 4 * vd1 ** 3
    ) * k ** 9

    a00 = 3.0 * wpi * wpi * (alpha0 * vthi0 * vthi0 + alpha1 * vthi1 * vthi1) * (vd0 * vd1) ** 4 * k ** 10

    return numpy.array(
        [a12, a11, a10, a09, a08, a07, a06, a05, a04, a03, a02, a01, a00],
        dtype=float,
    )


def wave_equation_two_electrons_twelveth_roots(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k
):
    k = numpy.asarray(k, dtype=float)
    w = numpy.empty((len(k), 12), dtype=complex)

    for i, ktmp in enumerate(k):
        coeff = wave_equation_two_electrons_twelveth_coefficients(
            wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, ktmp
        )
        w[i] = numpy.roots(coeff)

    return w


PAIR_CROSS = (0, 11)
PAIR_SELF = [(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]


def WE2Es12th_gfun(w, a=1.0, b=1.0):
    return a * numpy.real(w) + b * numpy.imag(w)


def WE2Es12th_track_branches_symmetric(
    k,
    w_all,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8,
):
    return _track_branches_symmetric(
        k,
        w_all,
        nm=12,
        cross_pair=PAIR_CROSS,
        self_pairs=PAIR_SELF,
        gfun=WE2Es12th_gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
        lam_sep=lam_sep,
        tol=tol,
        repair_passes=repair_passes,
    )


def wave_equation_two_electrons_twelveth_solve(
    wpe,
    mu,
    alpha0,
    alpha1,
    vd0,
    vd1,
    vthe0,
    vthe1,
    k,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8,
):
    w_raw = wave_equation_two_electrons_twelveth_roots(
        wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k
    )

    omega = WE2Es12th_track_branches_symmetric(
        k,
        w_raw,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
        lam_sep=lam_sep,
        tol=tol,
        repair_passes=repair_passes,
    )

    index=numpy.argsort(omega[:,-1])
    index=index[::-1]
    omega=omega[index,:]

    return omega


def WE2Es12th_check_symmetry(k, omega, a=1.0, b=1.0, tol=1.0e-12):
    return _check_symmetry(
        k,
        omega,
        cross_pair=PAIR_CROSS,
        self_pairs=PAIR_SELF,
        gfun=WE2Es12th_gfun,
        a=a,
        b=b,
        tol=tol,
    )


def WE2Es12th_check_local_slope_jumps(k, omega):
    return _check_local_slope_jumps(k, omega)


# ============================================================
# 6th-order coefficients
# optimized with scalar precomputation
# ============================================================
def two_stream_wave_equation_electron_sixth_coefficients(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k
):
    wpi = wpe / numpy.sqrt(mu)

    k2 = k * k
    k3 = k2 * k
    k4 = k2 * k2

    vd0_2 = vd0 * vd0
    vd1_2 = vd1 * vd1
    vd0vd1 = vd0 * vd1

    alpha_sum = alpha0 + alpha1
    wpi2 = wpi * wpi
    wpe2 = wpe * wpe

    a6 = 1.0

    a5 = -2.0 * k * (vd0 + vd1)

    a4 = k2 * (vd1_2 + 4.0 * vd0vd1 + vd0_2) - alpha_sum * (wpi2 + wpe2)

    a3 = (
        -2.0 * k3 * (vd0 * vd1_2 + vd0_2 * vd1)
        + 2.0 * k * (vd0 + vd1) * alpha_sum * wpi2
        + 2.0 * k * (alpha1 * vd0 + alpha0 * vd1) * wpe2
    )

    a2 = (
        k4 * vd0_2 * vd1_2
        - k2 * (vd1_2 + 4.0 * vd0vd1 + vd0_2) * alpha_sum * wpi2
        - k2 * (alpha0 * vd1_2 + alpha1 * vd0_2) * wpe2
    )

    a1 = 2.0 * k3 * (vd0 * vd1_2 + vd0_2 * vd1) * alpha_sum * wpi2

    a0 = -k4 * vd0_2 * vd1_2 * alpha_sum * wpi2

    return numpy.array([a6, a5, a4, a3, a2, a1, a0], dtype=float)


def two_stream_wave_equation_electron_sixth_roots(
    wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k
):
    k = numpy.asarray(k, dtype=float)
    w = numpy.empty((len(k), 6), dtype=complex)

    for i, ktmp in enumerate(k):
        coeff = two_stream_wave_equation_electron_sixth_coefficients(
            wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, ktmp
        )
        w[i] = numpy.roots(coeff)

    return w


PAIR_CROSS_6 = (0, 5)
PAIR_SELF_6 = [(1, 4), (2, 3)]


def WE2Es6th_gfun(w, a=1.0, b=1.0):
    return a * numpy.real(w) + b * numpy.imag(w)


def WE2Es6th_track_branches_symmetric(
    k,
    w_all,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8,
):
    return _track_branches_symmetric(
        k,
        w_all,
        nm=6,
        cross_pair=PAIR_CROSS_6,
        self_pairs=PAIR_SELF_6,
        gfun=WE2Es6th_gfun,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
        lam_sep=lam_sep,
        tol=tol,
        repair_passes=repair_passes,
    )


def two_stream_wave_equation_electron_sixth_solve(
    wpe,
    mu,
    alpha0,
    alpha1,
    vd0,
    vd1,
    vthe0,
    vthe1,
    k,
    a=1.0,
    b=1.0,
    lam_dist=1.0,
    lam_slope=1.0,
    lam_g=0.8,
    lam_gslope=0.8,
    lam_curv=0.25,
    lam_sep=0.15,
    tol=1.0e-12,
    repair_passes=8,
):
    w_raw = two_stream_wave_equation_electron_sixth_roots(
        wpe, mu, alpha0, alpha1, vd0, vd1, vthe0, vthe1, k
    )

    omega = WE2Es6th_track_branches_symmetric(
        k,
        w_raw,
        a=a,
        b=b,
        lam_dist=lam_dist,
        lam_slope=lam_slope,
        lam_g=lam_g,
        lam_gslope=lam_gslope,
        lam_curv=lam_curv,
        lam_sep=lam_sep,
        tol=tol,
        repair_passes=repair_passes,
    )


    index=numpy.argsort(omega[:,-1])
    index=index[::-1]
    omega=omega[index,:]

    return omega


def WE2Es6th_check_symmetry(k, omega, a=1.0, b=1.0, tol=1.0e-12):
    return _check_symmetry(
        k,
        omega,
        cross_pair=PAIR_CROSS_6,
        self_pairs=PAIR_SELF_6,
        gfun=WE2Es6th_gfun,
        a=a,
        b=b,
        tol=tol,
    )


def WE2Es6th_check_local_slope_jumps(k, omega):
    return _check_local_slope_jumps(k, omega)




############################################################
############################################################


def normal_pdf(x,mu,sigma):
    pdf= 1.0/numpy.sqrt(2.0*numpy.pi*sigma*sigma)*numpy.exp(-0.5*numpy.power((x-mu)/sigma,2.0))

    return pdf



def two_stream_growth_rate_electron_solve(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    #print(k,numpy.real(w),numpy.imag(w))

    wpi   = wpe/numpy.sqrt(mu)
    vthi0 = vthe0/numpy.sqrt(mu)
    vthi1 = vthe1/numpy.sqrt(mu)

    wr = numpy.real(w)
    H  = numpy.power(wpi/vthi0,2.0)*alpha0*wr*normal_pdf(wr/k,0.0,vthi0)+numpy.power(wpi/vthi1,2.0)*alpha1*wr*normal_pdf(wr/k,0.0,vthi1)+numpy.power(wpe/vthe0,2.0)*alpha0*(wr-k*vd0)*normal_pdf(wr/k,vd0,vthe0)+numpy.power(wpe/vthe1,2.0)*alpha1*(wr-k*vd1)*normal_pdf(wr/k,vd1,vthe1)
    H  = H/numpy.abs(k)

    #PG = 2.0*(alpha0+alpha1)*wpi*wpi*numpy.power(wr,-3.0)+12.0*numpy.power(k*wpi,2.0)*(alpha0*vthi0*vthi0+alpha1*vthi1*vthi1)*numpy.power(wr,-5.0)+2.0*alpha0*wpe*wpe*numpy.power(wr-k*vd0,-3.0)+12.0*numpy.power(k*wpe,2.0)*alpha0*vthe0*vthe0*numpy.power(wr-k*vd0,-5.0)+2.0*alpha1*wpe*wpe*numpy.power(wr-k*vd1,-3.0)+12.0*numpy.power(k*wpe,2.0)*alpha1*vthe1*vthe1*numpy.power(wr-k*vd1,-5.0)
    #PG = PG*numpy.power(k,2.0)

    PG1 = alpha0*numpy.power(wr,-3.0)+alpha1*numpy.power(wr,-3.0)+6.0*k*k*numpy.power(wr,-5.0)*(alpha0*vthi0*vthi0+alpha1*vthi1*vthi1)
    PG2 = alpha0*numpy.power(wr-k*vd0,-3.0)+alpha1*numpy.power(wr-k*vd1,-3.0)+6.0*k*k*alpha0*vthe0*vthe0*numpy.power(wr-k*vd0,-5.0)+6.0*k*k*alpha1*vthe1*vthe1*numpy.power(wr-k*vd1,-5.0)
    PG  = (PG1*wpi*wpi + PG2*wpe*wpe)*numpy.power(k,2.0)*2.0

    return -1.0*numpy.pi*H/PG



def two_stream_growth_rate_electron(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    gamma=[]
    for i in range(len(k)):
        gamma_tmp=two_stream_growth_rate_electron_solve(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k[i],w[i])
        gamma.append(gamma_tmp)
    
    return numpy.array(gamma)



############################################################
# electron-acoustic
############################################################



def wave_equation_EA(wpe,mu,alpha0,alpha1,vb,vthe0,vthe1,k):
    
    wpi  = wpe/numpy.sqrt(mu)
    vthi0=vthe0/numpy.sqrt(mu)
    vthi1=vthe1/numpy.sqrt(mu)

    Atmp = 1.0+alpha1*numpy.power(wpe/k/vthe1,2.0)
     
    Btmp = wpi*wpi+alpha0*wpe*wpe

    #Ctmp = 3.0*k*k*(alpha0*numpy.power(wpe*vthe0,2.0)+alpha0*numpy.power(wpi*vthi0,2.0)+alpha1*numpy.power(wpi*vthe1,2.0))

    Ctmp = 3.0*k*k*(alpha0*numpy.power(wpe*vthe0,2.0))


    xtmp=1.0+numpy.sqrt(1.0+4.0*Atmp*Ctmp/Btmp/Btmp)
    xtmp=0.5*Btmp/Atmp*xtmp

    return numpy.sqrt(xtmp)




def wave_equation_EA_gary(wpe,mu,alpha0,alpha1,vb,vthe0,vthe1,k):
    
    lambda0=vthe0/wpe
    lambda1=vthe1/wpe

    
    Atmp = alpha0*wpe*wpe*(1.0+3.0*k*k*lambda0*lambda0/alpha0)

    Btmp = 1.0+alpha1/lambda1/lambda1/k/k


    return numpy.sqrt(Atmp/Btmp)






############################################################
# ion-acoustic
############################################################

def wave_equation_two_electrons_IA_cs(wpe,mu,alpha0,alpha1,vthe0,vthe1):

    wpi     = wpe/numpy.sqrt(mu)

    vthe    = vthe0*vthe1/numpy.sqrt(alpha0*numpy.power(vthe1,2.0)+alpha1*numpy.power(vthe0,2.0))

    lambdaD = vthe/wpe

    cs      = wpi*lambdaD

    return cs


def wave_equation_two_electrons_IA_root(wpe,mu,alpha0,alpha1,vthe0,vthe1,k):

    wpi     = wpe/numpy.sqrt(mu)
    vthi0   = vthe0/numpy.sqrt(mu)
    vthi1   = vthe1/numpy.sqrt(mu)

    vthe    = vthe0*vthe1/numpy.sqrt(alpha0*numpy.power(vthe1,2.0)+alpha1*numpy.power(vthe0,2.0))

    lambdaD = vthe/wpe

    cs      = wpi*lambdaD

    wpow1    = 0.5*(alpha0+alpha1)*numpy.power(cs*k,2.0)/(1.0+numpy.power(k*lambdaD,2.0))
    wpow2    = (1.0+numpy.sqrt(1.0+12.0*(alpha0*vthi0*vthi0+alpha1*vthi1*vthi1)*(1.0+numpy.power(k*lambdaD,2.0))*numpy.power((alpha0+alpha1)*cs,-2.0)))

    return numpy.sqrt(wpow1*wpow2)


def wave_equation_two_electrons_IA_solve(wpe,mu,alpha0,alpha1,vthe0,vthe1,k):

    w=[]
    for ktmp in k:
        wtmp=wave_equation_two_electrons_IA_root(wpe,mu,alpha0,alpha1,vthe0,vthe1,ktmp)
        w.append(wtmp)
    
    return numpy.array(w)





def wave_equation_two_electrons_IA_growth_rate_root(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    #print(k,numpy.real(w),numpy.imag(w))

    wpi   = wpe/numpy.sqrt(mu)
    vthi0 = vthe0/numpy.sqrt(mu)
    vthi1 = vthe1/numpy.sqrt(mu)


    wr = numpy.real(w)
    H  = numpy.power(wpi/vthi0,2.0)*alpha0*wr*normal_pdf(wr/k,0.0,vthi0)+numpy.power(wpi/vthi1,2.0)*alpha1*wr*normal_pdf(wr/k,0.0,vthi1)+numpy.power(wpe/vthe0,2.0)*alpha0*(wr-k*vd0)*normal_pdf(wr/k,vd0,vthe0)+numpy.power(wpe/vthe1,2.0)*alpha1*(wr-k*vd1)*normal_pdf(wr/k,vd1,vthe1)

    H  = H/numpy.abs(k)

    #wr = numpy.real(w)
    #PG = 2.0*(alpha0+alpha1)*wpi*wpi*numpy.power(wr,-3.0)+12.0*numpy.power(k*wpi,2.0)*numpy.power(wr,-5.0)*(alpha0*vthi0*vthi0+alpha1*vthi1+vthi1)

    PG  = alpha0*numpy.power(wr,-3.0)*(1.0+6.0*numpy.power(k*vthi0/wr,2.0))+alpha1*numpy.power(wr,-3.0)*(1.0+6.0*numpy.power(k*vthi1/wr,2.0))
    PG  = PG*numpy.power(k*wpi,2.0)*2.0

    return -1.0*numpy.pi*H/PG


def wave_equation_two_electrons_IA_growth_rate_solve(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    gamma=[]
    for i in range(len(k)):
        gamma_tmp=wave_equation_two_electrons_IA_growth_rate_root(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k[i],w[i])
        gamma.append(gamma_tmp)
    
    return numpy.array(gamma)






def wave_equation_two_electrons_IA_root_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k):

    wpi     = wpe/numpy.sqrt(mu)

    vthe    = vthe0*vthe1/numpy.sqrt(alpha0*numpy.power(vthe1,2.0)+alpha1*numpy.power(vthe0,2.0))

    lambdaD = vthe/wpe

    cs      = wpi*lambdaD

    wpow    = (alpha0+alpha1)*numpy.power(cs*k,2.0)/(1.0+numpy.power(k*lambdaD,2.0))

    return numpy.sqrt(wpow)



def wave_equation_two_electrons_IA_solve_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,k):

    w=[]
    for ktmp in k:
        wtmp=wave_equation_two_electrons_IA_root_2(wpe,mu,alpha0,alpha1,vthe0,vthe1,ktmp)
        w.append(wtmp)
    
    return numpy.array(w)






def wave_equation_two_electrons_IA_growth_rate_root_2(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    #print(k,numpy.real(w),numpy.imag(w))

    wpi   = wpe/numpy.sqrt(mu)
    vthi0 = vthe0/numpy.sqrt(mu)
    vthi1 = vthe1/numpy.sqrt(mu)


    wr = numpy.real(w)
    H  = numpy.power(wpi/vthi0,2.0)*alpha0*wr*normal_pdf(wr/k,0.0,vthi0)+numpy.power(wpi/vthi1,2.0)*alpha1*wr*normal_pdf(wr/k,0.0,vthi1)+numpy.power(wpe/vthe0,2.0)*alpha0*(wr-k*vd0)*normal_pdf(wr/k,vd0,vthe0)+numpy.power(wpe/vthe1,2.0)*alpha1*(wr-k*vd1)*normal_pdf(wr/k,vd1,vthe1)

    H  = H/numpy.abs(k)

    #wr = numpy.real(w)
    PG = 2.0*(alpha0+alpha1)*wpi*wpi*numpy.power(wr,-3.0)
    PG = PG*numpy.power(k,2.0)

    return -1.0*numpy.pi*H/PG



def wave_equation_two_electrons_IA_growth_rate_solve_2(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k,w):

    gamma=[]
    for i in range(len(k)):
        gamma_tmp=wave_equation_two_electrons_IA_growth_rate_root_2(wpe,mu,alpha0,alpha1,vd0,vd1,vthe0,vthe1,k[i],w[i])
        gamma.append(gamma_tmp)
    
    return numpy.array(gamma)




############################################################
print("*"*65)
print("*"*65)