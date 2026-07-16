"""Python translation of the Christy-Bodek universal fit Fortran code that
calculates fit values for RL RT response functions and cross-section of
electron-nucleus scattering.

Fortran core source codes (07/06/26 archive):
    responseq.f
    csfitcomp.f
    gsmearing.f
    qenuc21off.f
    mec2021.f
    quasideut.f
    sf.f
    rescsp.f / rescsn.f / resmodp.f / resmodn.f
    nuc12sf.f / nuccs12cs.f / nucffs12c.f / nucffs12ct.f
    formfacts.f / vcoul.f

response output columns:
    qv, q2, ex, nu, rttot, rltot, rtqe, rlqe, rtie, rlie, rte, rle, rtns, rlns

cross-section output columns:
    e, theta, nu, ep, ex, q2, w2, q2_coulomb, w2_coulomb, epsilon_coulomb,
    flux_coulomb, x_coulomb, effective_potential, focusing_factor, xs_total,
    xs_qe, xs_inelastic, xs_mec, xs_narrow_states, xs_non_nuclear


Notes
-----
- use 1-based arrays for the fit parameter vectors (`xvalc`, resonance
  parameters, etc.) so the Python indexing matches the Fortran indexing.
- The algebra is kept close to the original code.
- Numerical guards are added only where the Fortran would otherwise use
  undefined/out-of-range values or encounter invalid square roots.
- The archive does not include the old cross-section driver program; the
  existing Python wrapper is retained, while VCOUL and NUCCS12CS are updated.
"""


from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------


def _pad1(values: Sequence[float]) -> List[float]:
    """Return a 1-based list: out[1] is the first physics parameter."""
    return [0.0] + list(values)


def _sqrt_pos(x: float) -> float:
    return math.sqrt(max(0.0, x))


def _clamp(x: float, lo: float, hi: float) -> float:
    return min(hi, max(lo, x))


def _exp(x: float) -> float:
    return math.exp(x)


# -----------------------------------------------------------------------------
# Constants and fit tables from the 07/06/26 Fortran source set
# -----------------------------------------------------------------------------

MP_MAIN = 0.938273
MP = 0.938272
MN = 0.939565
ALPHA = 1.0 / 137.036
PI = 3.141593
PI2 = PI * PI


XVALC = _pad1([
    0.18408, 6.9756, 0.17765, 12.545, 0.84099,
    1.212, 2.9829, 1.5953, 0.64928, -3.4309,
    58.794, 0.17722, 7.7061, 1, 0.21511,
    0.50422, 0.26811, 0.24094, 3.7362, 0.018371,
    0.080561, 0.15184, 0.222, 0, 0.222,
    0.046401, 0.66384, -0.035015, 0.022589, 0.086312,
    0.913, 0.13026, 0.16315, 1, 0.010472,
    6.9558, 2.1098, 0, 24.792, 0.01,
    1, 0.96116, 1.0726, 0.99999, 0.93939,
    1.0094, 0.96949, 1.0255, 0.98498, 0.99732,
    1.0414, 1, 1.0211, 1.035, 1.029,
    0.95411, 1.1187, 1.19, 1.1745, 1,
    1.0012, 1.0687, 1.0269, 1.1001, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
])


RESCSP_XVAL = _pad1([
    1.229, 1.5245, 1.5053, 1.713, 1.6809,
    1.4382, 0.12604, 0.235, 0.089515, 0.085287,
    0.075, 0.37589, 7.5995, 5.5341, 6.2953,
    0.58604, 6.7857, 2.4949, 0.43242, 1.7928,
    0.4976, 22.988, 0.020914, 0.92496, 0.029301,
    100.46, 0.22301, 0.49346, 2.3833, 5.266,
    14.292, 0.422, 0.05365, 104.36, 0.16219,
    0.55858, 0.27, 0.32575, 2.0218, 0.096124,
    0.71437, 1.1438, 0.50419, 1.6961, 0.36714,
    0.51125, 1.9923, 0.57, 4.2071, 0.71386,
    1, 0.99887, 0.99292, 1.0078, 0.9897,
    0.99921, 0.98542, 1.0247, 1.045, 0.99335,
    1.012, 0.99404, 7.7754, 3.2542, 1.8685,
    3.7473, 0, 93.168, 5.7939, 0.001021,
    3.1525, 7.3161, 2.2053, 0.53151, 3.2473,
    1.3766, 1.4511, 1.2058, 0.00021044, 0.79047,
    0.78882, 0.71543, 0.095268, 10.405, 1.8309,
    0, 13.763, 18.38, 0.056863, 0.016637,
    -0.20137, 5.9193, -0.19076, 0.40936, 547.81,
    111.13, 0, 0, 0, 0.07167,
])


RESCSN_XVAL = _pad1([
    1.229, 1.524, 1.5053, 1.7118, 1.6804,
    1.4387, 0.12659, 0.235, 0.089418, 0.086589,
    0.075, 0.37433, 6.9244, 3.31, 0.16188,
    1.0542, 4.809, 6.1513, 0.50963, 2.6077,
    0.33567, 64.341, 1.7105, 0.97942, 0.0029592,
    922.59, 0.29395, 0.74421, 0.0058413, 369.95,
    0.34929, 0.60553, 4.335, 89.525, 77726,
    7.5726, 0.25, 0.051415, 0.78543, 0.031374,
    0.13276, 0.50608, 0.57698, 2.614, 0.5539,
    0.34914, 1.93, 0.29, 2.2569, 0.24575,
    1.0138, 0.9846, 0.96512, 0.98146, 1.0425,
    1.0186, 0.96948, 0.98721, 1.0069, 1.0073,
    1.0091, 0.99948, 20.262, 0.15064, 3.513,
    0.88067, 27.222, 495.97, 26.483, 0.30933,
    1.6402e-08, 1.8731, 0.96153, 4.0622e-05, 2.5059e-09,
    2.4108, 1.3778, 0.01916, 4.9473, 7.5335,
    2.8901, 1.0402, 0.40764, 9.4392e-06, 0.073901,
    0.010134, 52.444, 3.8284, 0.0018283, 0.10426,
    -1.5743, 46.191, 51.493, 3.8696, 0.43902,
    0.97992, 0, 0, 1.0092, 0.96084,
])

# -----------------------------------------------------------------------------
# formfacts.f
# -----------------------------------------------------------------------------


def formfacts(q2: float) -> Tuple[float, float, float, float]:
    """FORMFACTS(q2, gmp, gep, gmn, gen)."""
    mu_p = 2.792782
    mu_n = -1.913148
    mp = 0.9382727
    tau = q2 / (4.0 * mp * mp)

    gmp = mu_p * (1.0 + 0.099481 * tau) / (
        1.0 + 11.089 * tau + 19.374 * tau * tau + 5.7798 * tau**3
    )
    gep = (1.0 + 0.24482 * tau * tau) / (
        1.0 + 11.715 * tau + 11.964 * tau * tau + 27.407 * tau**3
    )
    gd = (1.0 / (1.0 + q2 / 0.71)) ** 2

    gmn = mu_n * (1.0 + 2.330 * tau) / (
        1.0 + 14.720 * tau + 24.200 * tau**2 + 84.100 * tau**3
    )
    gen = (1.700 * tau / (1.0 + 3.300 * tau)) * gd

    gen = gen * ((q2 + 1189.4) / 1189.4) ** 219.73
    gmn = gmn / (((q2 + 0.35590) / 0.35590) ** 0.093020)
    return gmp, gep, gmn, gen


# -----------------------------------------------------------------------------
# resmodp.f / resmodn.f helpers (07/06/26 source set)
# -----------------------------------------------------------------------------


def _resmod_common(sf: int, w2: float, q2: float, xval: Sequence[float], *, proton: bool) -> float:
    """Shared implementation of the updated RESMODP and RESMODN routines."""
    if sf not in (1, 2):
        raise ValueError(f"sf must be 1 (T) or 2 (L), got {sf}")
    if w2 <= 0.0 or q2 < 0.0:
        return 0.0

    mp = 0.9382727 if proton else 0.939565
    mpi = 0.134977
    meta = 0.547862
    mp2 = mp * mp
    w = math.sqrt(w2)

    q20_seed = 0.05
    q20 = xval[50] * (1.0 + 0.25 * q20_seed / (q20_seed + q2)) ** 2.0
    log0 = math.log(q20 / 0.25**2)
    logq = math.log((q2 + q20) / 0.25**2)
    if log0 <= 0.0 or logq <= 0.0:
        return 0.0
    t = math.log(logq / log0)

    br = [[0.0] * 4 for _ in range(8)]
    ang = [0.0] * 8
    for i, v in enumerate([1.0, 0.45, 0.60, 0.65, 0.60, 0.65, 0.60], start=1):
        br[i][1] = v
    for i, v in enumerate([0.0, 0.40, 0.08, 0.0, 0.20, 0.0, 0.0], start=1):
        br[i][3] = v
    for i in range(1, 8):
        br[i][2] = 1.0 - br[i][1] - br[i][3]
    for i, v in enumerate([1.0, 0.0, 2.0, 3.0, 0.0, 1.0, 3.0], start=1):
        ang[i] = v

    x0 = [0.0] + [0.16] * 7
    if sf == 2:
        x0[1] = 0.015

    dip = 1.0 / (1.0 + q2 / 0.75) ** 2.0
    _mon = 1.0 / (1.0 + q2 / (1.0 if proton else 1.4))

    xb = q2 / (q2 + w2 - mp2)
    xpr1 = 1.0 / (1.0 + (w2 - (mp + mpi) ** 2) / (q2 + q20))
    xpr2 = 1.0 / (1.0 + (w2 - (mp + meta) ** 2) / (q2 + q20))
    if w <= mp + mpi:
        xpr1 = 1.0
    if w <= mp + meta:
        xpr2 = 1.0

    k = (w2 - mp2) / (2.0 * mp)
    kcm = (w2 - mp2) / (2.0 * w)
    if k == 0.0 or kcm == 0.0:
        return 0.0

    epicm = (w2 + mpi**2 - mp2) / (2.0 * w)
    ppicm = _sqrt_pos(epicm**2 - mpi**2)
    epi2cm = (w2 + (2.0 * mpi) ** 2 - mp2) / (2.0 * w)
    ppi2cm = _sqrt_pos(epi2cm**2 - (2.0 * mpi) ** 2)
    eetacm = (w2 + meta**2 - mp2) / (2.0 * w)
    petacm = _sqrt_pos(eetacm**2 - meta**2)

    mass = [0.0] * 8
    intwidth = [0.0] * 8
    width = [0.0] * 8
    num = 0
    for i in range(1, 7):
        num += 1
        mass[i] = xval[num]
    for i in range(1, 7):
        num += 1
        intwidth[i] = xval[num]
        width[i] = intwidth[i]
    mass[7] = xval[47]
    intwidth[7] = xval[48]
    width[7] = intwidth[7]

    kr = [0.0] * 8
    kcmr = [0.0] * 8
    pgam = [0.0] * 8
    pwid = [[0.0] * 4 for _ in range(8)]

    for i in range(1, 8):
        mi = mass[i]
        mi2 = mi * mi
        if mi <= 0.0:
            continue
        kr[i] = (mi2 - mp2) / (2.0 * mp)
        kcmr[i] = (mi2 - mp2) / (2.0 * mi)
        epicmr = (mi2 + mpi**2 - mp2) / (2.0 * mi)
        ppicmr = _sqrt_pos(epicmr**2 - mpi**2)
        epi2cmr = (mi2 + (2.0 * mpi) ** 2 - mp2) / (2.0 * mi)
        ppi2cmr = _sqrt_pos(epi2cmr**2 - (2.0 * mpi) ** 2)
        eetacmr = (mi2 + meta**2 - mp2) / (2.0 * mi)
        petacmr = _sqrt_pos(eetacmr**2 - meta**2)

        if ppicmr > 0.0:
            pwid[i][1] = intwidth[i] * (ppicm / ppicmr) ** (2.0 * ang[i] + 1.0) * (
                (ppicmr**2 + x0[i] ** 2) / (ppicm**2 + x0[i] ** 2)
            ) ** ang[i]
        if ppi2cmr > 0.0:
            pwid[i][2] = intwidth[i] * (ppi2cm / ppi2cmr) ** (2.0 * ang[i] + 4.0) * (
                (ppi2cmr**2 + x0[i] ** 2) / (ppi2cm**2 + x0[i] ** 2)
            ) ** (ang[i] + 2.0)
            pwid[i][2] *= w / mi
        if i in (2, 5) and petacmr > 0.0:
            pwid[i][3] = intwidth[i] * (petacm / petacmr) ** (2.0 * ang[i] + 1.0) * (
                (petacmr**2 + x0[i] ** 2) / (petacm**2 + x0[i] ** 2)
            ) ** ang[i]
        if kcmr[i] != 0.0:
            pgam[i] = intwidth[i] * (kcm / kcmr[i]) ** 2 * (
                (kcmr[i] ** 2 + x0[i] ** 2) / (kcm**2 + x0[i] ** 2)
            )
        width[i] = br[i][1] * pwid[i][1] + br[i][2] * pwid[i][2] + br[i][3] * pwid[i][3]

    height = [0.0] * 8
    rescoef = [[0.0] * 5 for _ in range(7)]
    for i in range(1, 7):
        for j in range(1, 5):
            num += 1
            rescoef[i][j] = xval[num]
        if sf == 1:
            h = rescoef[i][1] * (
                1.0 + rescoef[i][2] * q2 / (1.0 + rescoef[i][3] * q2)
            ) * dip ** rescoef[i][4]
            h /= (1.0 + 0.09 * q2) ** 2.0
        else:
            h = (rescoef[i][1] + rescoef[i][2] * q2) * _exp(-rescoef[i][3] * q2)
        height[i] = h * h

    if sf == 2:
        if proton:
            h7 = (xval[16] + xval[20] * q2) * _exp(-xval[24] * q2)
        else:
            h7 = (xval[44] + xval[45] * q2) * _exp(-xval[46] * q2)
    else:
        h7 = xval[49] / (1.0 + 0.03 * q2) * dip ** xval[45]
    height[7] = h7 * h7

    # Read but do not otherwise use the three groups of nonresonant coefficients,
    # matching the Fortran storage layout.
    for _i in range(1, 4):
        for _j in range(1, 5):
            num += 1
            _ = xval[num]

    sig_res = 0.0
    for i in range(1, 8):
        denom = (w2 - mass[i] ** 2) ** 2 + (mass[i] * width[i]) ** 2
        if denom <= 0.0 or intwidth[i] == 0.0:
            continue
        sigr = width[i] * pgam[i] / denom
        sigr = height[i] * kr[i] / k * kcmr[i] / kcm * sigr / intwidth[i]
        sig_res += sigr
    sig_res *= w
    if sf == 2:
        sig_res *= q2

    sig_nr = 0.0
    if sf == 1 and xpr1 < 1.0:
        a0 = xval[37] / (1.0 + q2 / xval[42])
        t1 = xval[38] + xval[39] * t ** xval[44]
        t2 = xval[40] + xval[41] * t ** xval[43]
        if xpr1 <= 1.0:
            sig_nr = 389.4 * a0 * (1.0 - xpr1) ** t1 * xpr1 ** t2
        if xpr2 <= 1.0:
            sig_nr += xval[46] * 389.4 * a0 * (1.0 - xpr2) ** t1 * xpr2 ** t2
    elif sf == 2 and xpr1 < 1.0:
        if proton:
            a0 = xval[37] / (1.0 + q2 / xval[39]) ** 2.0
            t1 = xval[38] / (1.0 + xval[40] * t ** xval[42])
            t2 = xval[41] / (1.0 + xval[43] * t) ** 2.0
        else:
            a0 = xval[37] / (1.0 + q2 / xval[39]) ** 2.0 * (1.0 + 85.0 * t) ** 2.0
            t1 = xval[38] / (1.0 + q2 / xval[40]) + xval[32] * math.log(q2 + xval[36])
            t2 = xval[41] / (1.0 + q2 / xval[42]) ** xval[43]
        if xpr1 <= 1.0:
            sig_nr = 389.4 * a0 * xb * (1.0 - xpr1) ** t1 * xpr1 ** t2

    sig = sig_res + sig_nr
    # The final threshold comparison in the Fortran is algebraically always false;
    # leave the computed result unchanged.
    return max(sig, 0.0)

# -----------------------------------------------------------------------------
# rescsp.f / rescsn.f / sf.f
# -----------------------------------------------------------------------------


def rescsp(w2: float, q2: float) -> Tuple[float, float]:
    """Return proton transverse and longitudinal reduced cross sections."""
    xval1 = [0.0] * 51
    xvalL = [0.0] * 51
    for i in range(1, 51):
        xval1[i] = RESCSP_XVAL[i]
        xvalL[i] = RESCSP_XVAL[50 + i]
        if i <= 12 or i in (47, 48):
            xvalL[i] = xval1[i]
    return (
        _resmod_common(1, w2, q2, xval1, proton=True),
        _resmod_common(2, w2, q2, xvalL, proton=True),
    )


def rescsn(w2: float, q2: float) -> Tuple[float, float]:
    """Return neutron transverse and longitudinal reduced cross sections."""
    xval1 = [0.0] * 51
    xvalL = [0.0] * 51
    for i in range(1, 51):
        xval1[i] = RESCSN_XVAL[i]
        xvalL[i] = RESCSN_XVAL[50 + i]
        if i <= 12 or i in (47, 48):
            xvalL[i] = xval1[i]
    return (
        _resmod_common(1, w2, q2, xval1, proton=False),
        _resmod_common(2, w2, q2, xvalL, proton=False),
    )


def sf(w2: float, q2: float) -> Tuple[float, float, float, float, float, float]:
    """Convert proton/neutron reduced cross sections to structure functions."""
    if q2 <= 0.0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    mp = 0.938272
    mp2 = mp * mp
    pi = 3.14159
    alpha = 1.0 / 137.03599
    x = q2 / (q2 + w2 - mp2)
    sigtp, siglp = rescsp(w2, q2)
    sigtn, sigln = rescsn(w2, q2)
    pref = abs(w2 - mp2) / (0.3894e3 * pi * pi * alpha * 8.0)
    f1p = sigtp * pref
    f1n = sigtn * pref
    flp = siglp * 2.0 * x * pref
    fln = sigln * 2.0 * x * pref
    denom = 1.0 + 4.0 * mp2 * x * x / q2
    f2p = (2.0 * x * f1p + flp) / denom
    f2n = (2.0 * x * f1n + fln) / denom
    return f1p, flp, f2p, f1n, fln, f2n

# -----------------------------------------------------------------------------
# gsmearing.f (Oct. 4, 2025 version in the 07/06/26 archive)
# -----------------------------------------------------------------------------


def _source_removal_energy(a: float, *, qe: bool = False) -> float:
    """Removal-energy constants explicitly defined by the current Fortran."""
    ia = int(round(a))
    if qe:
        if ia == 12:
            return 0.0165
        if ia in (27, 40):
            return 0.010
    else:
        if ia == 12:
            return 0.0165
        if ia in (27, 40):
            return 0.0085
    raise ValueError(f"The 07/06/26 Fortran only initializes removal energy for A=12, 27, or 40; got A={a}.")


def gsmearing(z: float, a: float, w2: float, q2: float, xvalc: Sequence[float]) -> Tuple[float, float, float]:
    """GSMEARING(Z,A,W2,Q2,...), translated from the 07/06/26 source."""
    if q2 <= 0.0:
        return 0.0, 0.0, 0.0
    nbins = 102
    nwid = 3.75
    bw = 2.0 * nwid / float(nbins)
    exmin = _source_removal_energy(a, qe=False)

    mp = 0.938272
    mp2 = mp * mp
    x = q2 / (q2 + w2 - mp2)
    nu = (w2 + q2 - mp2) / (2.0 * mp)
    if nu <= 0.0:
        return 0.0, 0.0, 0.0
    qv = math.sqrt(nu * nu + q2)
    kappa2 = 1.0 + 4.0 * mp2 * x * x / q2

    if a < 3.0:
        raise ValueError("The current GSMEARING source does not initialize kF/Es for A < 3.")
    kf = xvalc[25] + 0.00055 * (a / 12.0 * math.sqrt(a / 12.0) - 1.0)
    es = exmin

    norm = math.sqrt(PI) / 1.0027
    f1p = f1n = f2p = f2n = flp = fln = 0.0
    pf = 0.5 * kf
    pf2 = 1.5 * pf
    dw2dpf = 2.0 * qv
    dw2des = 2.0 * (nu + mp)

    for ism in range(1, nbins + 1):
        xxp = -nwid + bw * float(ism - 1)
        fyuse = bw / math.sqrt(2.0) / norm * _exp(-0.5 * xxp * xxp)
        fyuse *= 1.0 - xvalc[22] + xvalc[22] * xxp * xxp
        wsqp = w2 + xxp * pf * dw2dpf - es * dw2des
        wsqp2 = w2 + xxp * pf2 * dw2dpf - es * dw2des

        frac = 0.0
        # The source loop is DO j=1,1, so only WSQP is active.
        for j in range(1, 2):
            if j == 1:
                fract = 1.0 - frac
                w2t = wsqp
            else:
                fract = frac
                w2t = wsqp2
            if w2t <= 1.159:
                continue

            xt = q2 / (q2 + w2t - mp2)
            xp_denom = q2 + xvalc[36] / (1.0 + xvalc[39] * q2**2)
            xp = 1.0 / (1.0 + (w2t - mp2) / xp_denom)
            xt2 = max(xt, 0.0001)

            emcfac = (
                xvalc[21]
                - xvalc[11] * ((xt2 - xvalc[12]) ** 2.0) ** xvalc[13]
                + xvalc[14] * (1.0 - xp * xp) ** xvalc[15]
            )
            emcfac /= (1.0 + 0.05 * xp) ** xvalc[40]
            emcfac_l = 1.0

            f1pp, flpp, f2pp, f1nn, flnn, f2nn = sf(w2t, q2)
            f1pp *= emcfac
            f1nn *= emcfac
            flpp *= emcfac * emcfac_l
            flnn *= emcfac * emcfac_l
            denom = 1.0 + 4.0 * xt * xt * mp2 / q2
            f2pp = (2.0 * xt * f1pp + flpp) / denom
            f2nn = (2.0 * xt * f1nn + flnn) / denom

            f1p += f1pp * fyuse * fract
            f1n += f1nn * fyuse * fract
            f2p += f2pp * fyuse * fract
            f2n += f2nn * fyuse * fract
            flp += flpp * fyuse * fract
            fln += flnn * fyuse * fract

    f1 = z * f1p + (a - z) * f1n
    f2 = z * f2p + (a - z) * f2n
    fl = z * flp + (a - z) * fln

    cof = _clamp(1.0 - (es / nu) ** 4, 0.0, 1.0)
    if nu <= es:
        cof = 0.0
    f1 *= cof
    f2 *= cof
    fl *= cof
    f2 = (2.0 * x * f1 + fl) / kappa2
    return max(f1, 0.0), max(f2, 0.0), max(fl, 0.0)

# -----------------------------------------------------------------------------
# qenuc21off.f
# -----------------------------------------------------------------------------


def qenuc21off(z: float, a: float, q2: float, w2: float, xvalc: Sequence[float]) -> Tuple[float, float]:
    """Updated superscaling quasielastic F1 and F2 per nucleus."""
    mp = 0.938272
    psimax = 5.0
    moff = mp
    ia = int(a)
    if ia == 1 or q2 <= 0.0:
        return 0.0, 0.0

    nu = (w2 - mp**2 + q2) / (2.0 * mp)
    if nu <= 0.0:
        return 0.0, 0.0
    tau = q2 / (4.0 * mp**2)
    qv = math.sqrt(nu * nu + q2)
    nuel = q2 / (2.0 * 0.931494 * a)
    ex = nu - nuel
    _ = ex

    gmp, gep, gmn, gen = formfacts(q2)
    if ia < 3:
        raise ValueError("The current QENUC21OFF source does not initialize kF/Es for A < 3.")
    esmin = _source_removal_energy(a, qe=True)
    kf = xvalc[23]
    es = esmin
    qvmax = 1.5
    if qv < qvmax:
        es += xvalc[24]
        es -= xvalc[24] * math.sqrt(max(0.0, 1.0 - (qv / qvmax) ** 2))
    nup = nu - es

    if qv > 2.0 * kf:
        pauli_sup2 = 1.0
    else:
        pauli_sup2 = 0.75 * (qv / kf) * (1.0 - (qv / kf) ** 2 / 12.0)
    _ = pauli_sup2

    kappa = qv / (2.0 * mp)
    lam = nu / (2.0 * mp)
    lamp = nup / (2.0 * mp)
    lampn = -lamp
    taup = kappa**2 - lamp**2
    xi = math.sqrt(1.0 + (kf / moff) ** 2) - 1.0
    if taup <= 0.0 or xi <= 0.0 or 1.0 + taup <= 0.0:
        return 0.0, 0.0

    root = math.sqrt(taup * (1.0 + taup))
    args = [
        (1.0 + lam) * taup + kappa * root,
        (1.0 + lamp) * taup + kappa * root,
        (1.0 + lampn) * taup + kappa * root,
    ]
    if min(args) <= 0.0:
        return 0.0, 0.0
    sqrt_xi = math.sqrt(xi)
    psi = (lam - taup) / sqrt_xi / math.sqrt(args[0])
    psip = (lamp - taup) / sqrt_xi / math.sqrt(args[1])
    psipn = (lampn - taup) / sqrt_xi / math.sqrt(args[2])

    nul = (q2 / qv / qv) ** 2
    nut = tau / (2.0 * kappa**2)
    avgn = a - z
    gm2bar = z * gmp**2 + avgn * gmn**2
    ge2bar = z * gep**2 + avgn * gen**2
    f1ff = (tau * math.sqrt(gm2bar) + math.sqrt(ge2bar)) / (1.0 + tau)
    f2ff = (math.sqrt(gm2bar) - math.sqrt(ge2bar)) / (1.0 + tau)
    geoff = f1ff - tau * moff / mp * f2ff
    gmoff = f1ff + moff / mp * f2ff

    delta = tau / kappa / kappa * xi * (1.0 - psi**2) * (
        kappa * math.sqrt(1.0 + 1.0 / tau) + xi / 3.0 * (1.0 - psi**2)
    )
    gl = kappa**2 / tau * (geoff**2 + (geoff**2 + tau * gmoff**2) * delta / (1.0 + tau))
    gt = 2.0 * tau * gmoff**2 + (geoff**2 + tau * gmoff**2) * delta / (1.0 + tau)
    norm = 2.0 * kappa * (1.0 + xi * (1.0 + psi**2) / 2.0)
    gl /= norm
    gt /= norm

    def fitted_scaling(p: float) -> float:
        value = xvalc[7] / (1.0 + xvalc[8] ** 2 * (p + xvalc[9]) ** 2)
        value /= 1.0 + _exp(xvalc[10] * p)
        value *= (1.0 - abs(p) / psimax) ** 2.0 * (1.0 + abs(p) / psimax) ** 2.0
        value /= 1.0 + 0.125 * p * p
        if p > psimax:
            return 0.0
        return max(0.0, value)

    fyp = fitted_scaling(psip)
    fyn = fitted_scaling(psipn)
    fy = max(0.0, fyp - fyn)

    x2 = qv / kf
    d1 = 1.0 + xvalc[31] * math.sqrt(a / 12.0)
    g1 = d1 * xvalc[28]
    g2 = xvalc[29]
    g3 = xvalc[30]
    if x2 > 4.0:
        pb2l = 1.0
    else:
        pb2l = 1.0 - g1 * (4.0 - x2) ** 2.5 - g2 * (4.0 - x2) ** 3.5
        pb2l -= g3 * (4.0 - x2) ** 1.5
        denom = (x2 - 0.18) ** 2.0
        pb2l = pb2l * (x2 - 0.2) ** 2 / denom if denom != 0.0 else 0.0
        pb2l = _clamp(pb2l, 0.0, 1.0)
    if psip > psimax:
        fy = 0.0
    fyl = pb2l * fy

    f2 = nu / kf * (fyl * nul * gl + fy * nut * gt)
    f1 = mp * fy / kf * gt / 2.0
    cof = _clamp(1.0 - (es / nu) ** 4.0, 0.0, 1.0)
    if nu <= es:
        cof = 0.0
    return max(f1 * cof, 0.0), max(f2 * cof, 0.0)

# -----------------------------------------------------------------------------
# mec2021.f / quasideut.f
# -----------------------------------------------------------------------------


def mec2021(z: float, a: float, w2: float, q2: float, xvalm: Sequence[float]) -> float:
    """Updated transverse-enhancement/MEC contribution F1."""
    mp = 0.938272
    mp2 = mp * mp
    if w2 <= 0.0 or q2 < 0.0 or a < 2.5:
        return 0.0
    nu = (w2 - mp2 + q2) / (2.0 * mp)
    if nu <= 0.0:
        return 0.0
    numin = _source_removal_energy(a, qe=False)
    nuel = q2 / (2.0 * 0.931494 * a)
    ex = nu - nuel
    q20 = xvalm[26]

    a1 = xvalm[1]
    a2 = xvalm[35]
    y = (q2 + q20) ** xvalm[2] / (xvalm[3] + q2) ** xvalm[4]
    y2 = y * (q2 + xvalm[18]) ** xvalm[19] / (xvalm[20] + q2) ** 2.0
    b1 = xvalm[5]
    b2 = xvalm[6]
    c1 = xvalm[32] + xvalm[33] * q2 * math.sqrt(q2)
    c2 = xvalm[17]
    if c1 == 0.0 or c2 == 0.0:
        return 0.0

    w2min = mp2 + 2.0 * mp * numin - q2
    w2min2 = 1.059
    dw2 = max(w2 - w2min, 0.0)
    dw22 = max(w2 - w2min2, 0.0)
    t1 = (w2 - b1) ** 2 / (2.0 * c1**2)
    f1mec = y * _exp(-t1) * dw2
    f1mec2 = y2 / ((w2 - b2) ** 2 + b2**2 * c2**2) * dw22
    f1mec = a * (a1 * f1mec + a2 * f1mec2)

    if nu < numin:
        f1mec = 0.0
    cof = _clamp(1.0 - (numin / nu) ** 4, 0.0, 1.0)
    if ex <= numin:
        cof = 0.0
    f1mec *= cof
    return f1mec if f1mec > 1.0e-9 else 0.0


def quasideut(z: float, a: float, w2: float, q2: float, xvalm: Sequence[float]) -> float:
    """Translate QUASIDEUT. The active 07/06/26 source explicitly sets sigqd=0."""
    _ = (z, a, w2, q2, xvalm)
    return 0.0

# -----------------------------------------------------------------------------
# nucffs12c.f / nucffs12ct.f / nuc12sf.f
# -----------------------------------------------------------------------------


def nucffs12c(a: float, z: float, q32: float, state: int) -> float:
    """Longitudinal 12C nuclear form factors from the updated source."""
    if q32 < 0.0:
        return 0.0
    q2f = q32 / (0.1975 * 0.1975)
    qf = math.sqrt(q2f)
    ff2 = 0.0
    fflowq2 = q2f**3.0 / (q2f**3.0 + 0.01) if q2f > 0.0 else 0.0

    if state == 1:
        # Q2 is undeclared in the Fortran state-1 branch. q32 is the intended
        # momentum-transfer argument and makes this otherwise-unused branch deterministic.
        radius = 2.45
        x2 = (q32 / 0.197328**2) * radius**2
        alp = (z - 2.0) / 3.0
        char = x2 * (2.0 + 3.0 * alp) / (12.0 + 30.0 * alp)
        ff = _exp(-char) * (1.0 - alp * x2 / (6.0 + 15.0 * alp)) if char < 80.0 else 0.0
        _ = ff
        alph = 4.0 / 3.0
        a0 = 1.65
        g0 = 8.0e-5 * _exp(-((q2f - 2.9) / 0.44) ** 2.0)
        h0 = (1.0 - alph * q2f * a0 * a0 / 2.0 / (2.0 + 3.0 * alph)) * _exp(-q2f * a0 * a0 / 4.0)
        g4 = 1.0e-5 * _exp(-((q2f - 4.0) / 1.2) ** 2.0)
        ff2 = h0 * h0 + g0 + g4
    elif state == 2:
        g1 = 1.41e-2 * _exp(-(q2f - 1.125) ** 2 / 1.71**2)
        g2 = 7.0e-4 * _exp(-(q2f - 3.7) ** 2 / 1.6**2)
        g3 = 3.3e-6 * _exp(-(q2f - 6.5) ** 2 / 7.0**2)
        ff2 = fflowq2 * (g1 + g2 + g3)
    elif state == 3:
        b = 1.3457
        terms = (
            0.52 * (b * qf) ** 2
            - 0.025 * (b * qf) ** 4
            - 0.7e-2 * (b * qf) ** 6
            + 0.5e-3 * (b * qf) ** 8
            - 0.5e-4 * (b * qf) ** 10
        )
        ff2 = (1.0 / z) * _exp(-0.5 * b * b * q2f) * terms if z != 0.0 else 0.0
        ff2 *= ff2
    elif state == 4:
        g1 = 5.0e-3 * _exp(-(q2f - 1.46) ** 2 / 1.6**2)
        g2 = 6.6e-4 * _exp(-(q2f - 3.46) ** 2 / 2.0**2)
        g3 = 7.0e-6 * _exp(-(q2f - 7.0) ** 2 / 2.8**2)
        ff2 = q2f**3 / (q2f**3 + 0.2) * (g1 + g2 + g3) if q2f > 0.0 else 0.0
    elif state == 5:
        ff2 = (5.0e-4 * _exp(-(qf - 1.0) ** 2 / 0.3**2) + 8.0e-4 * _exp(-(qf - 1.4) ** 2 / 0.4**2)) * fflowq2
    elif state == 6:
        g1 = 4.0e-4 * _exp(-(qf - 1.0) ** 2 / 0.35**2)
        g2 = 8.0e-4 * _exp(-(qf - 1.75) ** 2 / 0.45**2)
        g3 = 4.0e-4 * _exp(-(qf - 0.85) ** 2 / 0.65**2)
        # g4 is used without an active assignment in the Fortran; the commented
        # assignment and surrounding pattern indicate the intended value is zero.
        ff2 = (g1 + g2 + g3) * fflowq2
    elif state == 7:
        ff2 = 6.0e-4 * _exp(-(qf - 0.85) ** 2 / 0.7**2)
    elif state == 8:
        ff2 = 12.0e-4 * _exp(-(qf - 1.05) ** 2 / 0.6**2) * fflowq2
    elif state == 9:
        ff2 = 3.2e-4 * _exp(-(qf - 1.3) ** 2 / 0.5**2) * fflowq2
    elif state == 10:
        ff2 = (1.6e-4 * _exp(-(qf - 1.2) ** 2 / 0.42**2) + 1.6e-5 * _exp(-(qf - 1.8) ** 2 / 0.4**2)) * fflowq2
    elif state == 11:
        ff2 = (2.8e-3 * _exp(-(qf - 0.60) ** 2 / 0.15**2) + 6.9e-3 * _exp(-(qf - 0.84) ** 2 / 0.55**2)) * fflowq2
    elif state == 12:
        ff2 = 0.0047 * _exp(-(qf - 1.0) ** 2 / 0.48**2) * fflowq2
    elif state == 13:
        ff2 = 0.0029 * _exp(-(qf - 1.69) ** 2 / 0.7**2) * fflowq2

    if qf * qf > 12.0:
        ff2 = 0.0
    return math.sqrt(max(ff2, 0.0))


def nucffs12ct(a: float, z: float, q2: float, state: int) -> float:
    """Transverse 12C nuclear form factors from the updated source."""
    _ = (a, z)
    if q2 < 0.0:
        return 0.0
    q2f = q2 / (0.1975 * 0.1975)
    qf = math.sqrt(q2f)
    ff2 = 0.0
    if state == 14:
        ff2 = 2.5e-4 * _exp(-(qf - 0.63) ** 2 / 0.4**2) + 2.8e-4 * _exp(-(qf - 0.84) ** 2 / 0.2**2) - 2.5e-5 * _exp(-qf)
    elif state == 15:
        ff2 = 5.9e-4 * _exp(-(qf - 1.2) ** 2 / 0.55**2) + 2.4e-4 * _exp(-(qf - 2.2) ** 2 / 0.6**2)
    elif state == 16:
        ff2 = 2.6e-4 * _exp(-(qf - 1.6) ** 2 / 0.6**2) + 5.0e-5 * _exp(-(qf - 2.5) ** 2 / 0.35**2)
    elif state == 17:
        ff2 = 2.1e-4 * _exp(-(qf - 0.8) ** 2 / 0.35**2) + 1.6e-4 * _exp(-(qf - 1.2) ** 2 / 0.425**2)
    elif state == 18:
        ff2 = (
            9.5e-4 * _exp(-(qf - 1.27) ** 2 / 0.77**2)
            + 3.5e-4 * _exp(-(qf - 1.7) ** 2 / 0.6**2)
            + 1.0e-4 * _exp(-(qf - 2.2) ** 2 / 0.3**2)
            - 3.6e-4 * _exp(-qf)
        )
    elif state == 19:
        ff2 = 2.1e-4 * _exp(-(qf - 1.45) ** 2 / 0.5**2) + 5.5e-5 * _exp(-(qf - 2.1) ** 2 / 0.4**2)
    elif state == 20:
        ff2 = 1.8e-3 * _exp(-(qf - 0.8) ** 2 / 0.36**2) + 1.0e-4 * _exp(-(qf - 1.5) ** 2 / 0.5**2)
    elif state == 21:
        # g2 is not assigned in the source; the commented expression and active
        # pattern imply zero.
        ff2 = 9.0e-4 * _exp(-(qf - 0.35) ** 2 / 0.3**2)
    if qf * qf > 12.0:
        ff2 = 0.0
    return math.sqrt(max(ff2, 0.0))


def nuc12sf(z: float, a: float, nu: float, q2p: float, state: int) -> Tuple[float, float]:
    """Narrow-state F1 and FL used by RESPONSEQ."""
    if not (1 <= state <= 21) or nu <= 0.0:
        # RESPONSEQ loops through 22 although the source arrays contain 21 states.
        return 0.0, 0.0
    mp = 0.93827
    exc = [
        0.0,
        0.0, 0.00444, 0.00765, 0.00964, 0.01084, 0.0137, 0.0151,
        0.0161, 0.0183, 0.020, 0.0230, 0.0315, 0.042, 0.0151,
        0.0161, 0.0166, 0.0181, 0.0193, 0.0206, 0.0235, 0.0315,
    ]
    wid = [
        0.0,
        0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00125,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00475, 0.009, 0.012,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00035, 0.00015, 0.004, 0.009,
    ]
    q2 = max(q2p, 0.0)
    qv2 = q2 + nu * nu
    # smwid = 0.0035
    smwid = 0.001

    width = math.sqrt(smwid * smwid + wid[state] * wid[state])
    norm = width * math.sqrt(3.14159)
    nuel = q2 / (2.0 * 0.931494 * a)
    q2f = q2 / (0.1975 * 0.1975)
    if state == 18:
        exc[state] = min(0.0194 + 0.00016 * math.sqrt(q2f), 0.01955)
    nuex = nuel + exc[state]
    fs = _exp(-(nu - nuex) ** 2 / width**2) / norm
    if (nu - nuex) / width > 5.0:
        fs = 0.0

    ff = nucffs12c(a, z, qv2, state) if state <= 13 else 0.0
    fft = nucffs12ct(a, z, qv2, state) if state > 13 else 0.0
    w1 = 0.5 * (z * fft) ** 2 * fs
    wl = (z * ff) ** 2 * fs
    f1 = mp * w1
    fl = q2 * q2 / qv2 / nu * wl if qv2 > 0.0 else 0.0
    return f1, fl

# -----------------------------------------------------------------------------
# csfitcomp.f
# -----------------------------------------------------------------------------


def csfitcomp(w2: float, q2: float, a: float, z: float, xvalc: Sequence[float], kind: int) -> Tuple[float, float]:
    """Return (sigma_T, sigma_L) for the requested component type."""
    if q2 <= 0.0:
        return 0.0, 0.0
    psimin = -2.3
    psimax = 5.0
    nbins = 220
    mp = 0.938272
    mp2 = mp * mp
    x = q2 / abs(w2 - mp2 + q2)
    if x == 0.0:
        return 0.0, 0.0
    dpsi = (psimax - psimin) / float(nbins)
    kappa2 = 1.0 + 4.0 * x * x * mp2 / q2

    int2 = 0.0
    for i in range(1, nbins + 1):
        psip = psimin + dpsi * (i - 1)
        fy2 = xvalc[7] / (1.0 + xvalc[8] ** 2 * (psip + xvalc[9]) ** 2)
        fy2 /= 1.0 + _exp(xvalc[10] * psip)
        fy2 *= (1.0 - abs(psip) / psimax) ** 2.0 * (1.0 + abs(psip) / psimax) ** 2.0
        fy2 /= 1.0 + 0.08 * psip**2
        if psip > psimax:
            fy2 = 0.0
        int2 += max(0.0, fy2)
    if int2 == 0.0:
        return 0.0, 0.0
    rat = 1.0 / int2 / dpsi

    f1i, f2i, fli = gsmearing(z, a, w2, q2, xvalc)
    fli = max(fli, 0.0)
    f1qe, f2qe = qenuc21off(z, a, q2, w2, xvalc)
    f1qe *= rat
    f2qe *= rat
    flqe = max(kappa2 * f2qe - 2.0 * x * f1qe, 0.0)

    f1mec = mec2021(z, a, w2, q2, xvalc)
    flmec = 0.0
    f2mec = 2.0 * x * f1mec / kappa2
    f1qd = quasideut(z, a, w2, q2, xvalc)
    flqd = 0.0
    f2qd = 2.0 * x * f1qd / kappa2

    if kind == 1:
        f1 = f1i + f1qe + f1mec + f1qd
        f2 = f2i + f2qe + f2mec + f2qd
        fl = fli + flqe
    elif kind == 2:
        f1, f2, fl = f1qe, f2qe, flqe
    elif kind == 3:
        f1, f2, fl = f1i, f2i, fli
    elif kind == 4:
        f1, f2, fl = f1mec, f2mec, flmec
    elif kind == 5:
        f1, f2, fl = f1qe + f1mec, f2qe + f2mec, flqe + flmec
    else:
        raise ValueError(f"Unsupported CSFITCOMP type={kind}")
    _ = f2
    return f1, max(fl, 0.0) / (2.0 * x)

# -----------------------------------------------------------------------------
# Minimal RL/RT evaluation API
# -----------------------------------------------------------------------------

OUTPUT_COLUMNS = [
    "qv", "q2", "ex", "nu",
    "rttot", "rltot", "rtqe", "rlqe",
    "rtie", "rlie", "rte", "rle",
    "rtns", "rlns"
]


def calculate_response_point(qv: float, nu: float, *, a: float = 12.0, z: float = 6.0, xvalc: Sequence[float] = XVALC) -> Optional[dict]:
    """Calculate the updated RESPONSEQ response components at one point."""
    q2 = qv * qv - nu * nu
    if q2 <= 0.0 or nu <= 0.0:
        return None
    nuel = q2 / (2.0 * 0.931494 * a)
    ex = nu - nuel
    w2 = MP_MAIN * MP_MAIN + 2.0 * MP_MAIN * nu - q2
    xb = q2 / (2.0 * MP_MAIN * nu)

    def response(kind: int) -> tuple[float, float]:
        f1, sigl = csfitcomp(w2, q2, a, z, xvalc, kind)
        fl = 2.0 * xb * sigl
        rt = 2.0 / MP_MAIN * f1 / 1000.0
        rl = qv * qv / q2 / 2.0 / MP_MAIN / xb * fl / 1000.0
        return rt, rl

    rttot, rltot = response(1)
    rtqe, rlqe = response(2)
    rtie, rlie = response(3)
    rte, _rle_model = response(4)
    rle = 0.0

    f1ns = 0.0
    flns = 0.0
    for state in range(2, 23):
        f1_state, fl_state = nuc12sf(z, a, nu, q2, state)
        f1ns += f1_state
        flns += fl_state
    rtns = 2.0 / MP_MAIN * f1ns / 1000.0
    rlns = qv * qv / q2 / 2.0 / MP_MAIN / xb * flns / 1000.0
    rttot += rtns
    rltot += rlns

    return {
        "qv": qv, "q2": q2, "ex": ex, "nu": nu,
        "rttot": rttot, "rltot": rltot,
        "rtqe": rtqe, "rlqe": rlqe,
        "rtie": rtie, "rlie": rlie,
        "rte": rte, "rle": rle,
        "rtns": rtns, "rlns": rlns,
    }


def calculate_response_table(table: pd.DataFrame | np.ndarray | Iterable[tuple[float, float]], *, a: float = 12.0, z: float = 6.0, xvalc: Sequence[float] = XVALC) -> pd.DataFrame:
    """Calculate RL/RT responses for many kinematic points.

    Parameters
    ----------
    table : DataFrame, ndarray, or iterable of (qv, nu)
        Input kinematic table. If a DataFrame is passed, it must contain the
        columns ``qv`` and ``nu``.
    a : float
        Nuclear mass number A.
    z : float
        Nuclear charge Z.
    xvalc : sequence of float
        Fit parameter table used by the original Fortran code.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns in ``OUTPUT_COLUMNS`` order.
    """
    if isinstance(table, pd.DataFrame):
        qv_values = table["qv"].to_numpy(dtype=float)
        nu_values = table["nu"].to_numpy(dtype=float)
    else:
        arr = np.asarray(list(table) if not isinstance(table, np.ndarray) else table, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("Input table must have shape (N, 2) with columns [qv, nu].")
        qv_values = arr[:, 0]
        nu_values = arr[:, 1]

    rows = []
    for qv, nu in zip(qv_values, nu_values):
        row = calculate_response_point(float(qv), float(nu), a=a, z=z, xvalc=xvalc)
        if row is not None:
            rows.append(row)
        else:
            rows.append(
                {"qv": qv,
                "q2": 0.0,
                "ex": 0.0,
                "nu": nu,
                "rttot": 0.0,
                "rltot": 0.0,
                "rtqe": 0.0,
                "rlqe": 0.0,
                "rtie": 0.0,
                "rlie": 0.0,
                "rte": 0.0,
                "rle": 0.0,
                "rtns": 0.0,
                "rlns": 0.0})

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


# -----------------------------------------------------------------------------
# Cross-section driver logic from qemodplot_norm.f / vcoul.f / nuccs12cs.f
# -----------------------------------------------------------------------------


def vcoul(a: float, z: float) -> float:
    """Average Coulomb potential from VCOUL(A,Z,V), in GeV."""
    hbarc = 0.197327
    alpha = 1.0 / 137.0
    c_aste = 0.775
    r0 = 1.1 * a ** (1.0 / 3.0) + 0.86 * a ** (-1.0 / 3.0)
    v0 = 1.5 * alpha * hbarc * (z - 1.0) / r0
    return c_aste * v0


def nuccs12cs(z: float, a: float, beam_energy: float, scattered_energy: float,
              scattering_angle_deg: float, state: int) -> float:
    """Narrow-state cross section in microbarn, updated to smwid=3.5 MeV."""
    if not (1 <= state <= 22) or beam_energy <= 0.0:
        return 0.0
    mp = 0.93827
    alpha = 7.29735e-3
    pi = 3.14159
    radcon = 0.0174533
    exc = [
        0.0, 0.0, 0.00444, 0.00765, 0.00964, 0.01084, 0.0137, 0.0151,
        0.0161, 0.0183, 0.020, 0.0230, 0.0315, 0.042, 0.0151, 0.0161,
        0.0166, 0.0181, 0.0193, 0.0206, 0.0235, 0.0315, 0.01271,
    ]
    wid = [
        0.0, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00125,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00475, 0.009, 0.012,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00035, 0.00015, 0.004,
        0.009, 0.0002,
    ]
    # smwid = 0.0035
    smwid = 0.001

    width = math.sqrt(smwid * smwid + wid[state] * wid[state])
    norm = width * math.sqrt(pi)

    sin_half = math.sin(radcon * scattering_angle_deg / 2.0)
    sin2 = sin_half * sin_half
    cos2 = 1.0 - sin2
    if sin2 <= 0.0 or cos2 <= 0.0:
        return 0.0
    tan2 = sin2 / cos2
    nu = beam_energy - scattered_energy
    q2 = 4.0 * beam_energy * scattered_energy * sin2
    q2v = q2 + nu * nu
    if q2 <= 0.0 or q2v <= 0.0:
        return 0.0

    epnuc = a * 0.931494 * beam_energy / (a * 0.931494 + 2.0 * beam_energy * sin2)
    nuel = beam_energy - epnuc
    q2f = q2 / (0.1975 * 0.1975)
    state_exc = exc[state]
    if state == 18:
        state_exc = min(0.0193 + 0.00015 * math.sqrt(q2f), 0.0195)
    nuex = nuel + state_exc
    q2s = 4.0 * beam_energy * (beam_energy - nuex) * sin2
    q2vs = q2s + nuex * nuex
    if q2vs < 0.0:
        return 0.0

    fs = _exp(-(nu - nuex) ** 2 / width**2) / norm
    if (nu - nuex) / width > 4.0:
        fs = 0.0
    ff = nucffs12c(a, z, q2vs, state) if state <= 13 else 0.0
    fft = nucffs12ct(a, z, q2vs, state) if state > 13 else 0.0
    wl2 = (z * ff) ** 2
    wt2 = (z * fft) ** 2

    mott = 0.3894e3 * alpha**2 / 4.0 * cos2 / beam_energy**2 / sin2**2
    recoil = a * mp / 1.007276 / (12.0 * mp + beam_energy * (1.0 - math.cos(radcon * scattering_angle_deg)))
    sigma = 1000.0 * mott * recoil * (
        q2 * q2 / q2v / q2v * wl2 + (q2 / 2.0 / q2v + tan2) * wt2
    )
    return fs * sigma


def calculate_cross_section_point(
    beam_energy: float,
    scattering_angle_deg: float,
    energy_transfer: float,
    *,
    a: float = 12.0,
    z: float = 6.0,
    xvalc: Sequence[float] = XVALC,
    coulomb_correction: bool = True,
) -> Optional[dict]:
    """Calculate the inclusive cross section and its component pieces.

    Parameters
    ----------
    beam_energy : float
        Incoming electron energy E in GeV.
    scattering_angle_deg : float
        Electron scattering angle in degrees.
    energy_transfer : float
        Energy transfer nu in GeV.
    a : float, optional
        Nuclear mass number A. Default is 12 for carbon.
    z : float, optional
        Nuclear charge Z. Default is 6 for carbon.
    xvalc : sequence of float, optional
        1-based fit-parameter table used by `csfitcomp`.
    coulomb_correction : bool, optional
        If True, apply the same Coulomb energy shift and focusing factor used
        in `qemodplot_norm.f`.

    Returns
    -------
    dict or None
        A dictionary of kinematics and cross-section pieces, or None when the
        point is outside the region printed by the original Fortran driver.
    """
    mp = 0.9382727
    mp2 = mp * mp
    alpha = 1.0 / 137.0
    pi = 3.14159
    pi2 = pi * pi
    radcon = 0.0174533

    scattered_energy = beam_energy - energy_transfer
    if scattered_energy <= 0.0:
        return None

    sin_half_theta = math.sin(radcon * scattering_angle_deg / 2.0)
    sin2 = sin_half_theta * sin_half_theta
    cos2 = 1.0 - sin2
    if cos2 <= 0.0 or sin2 <= 0.0:
        return None
    tan2 = sin2 / cos2

    q2 = 4.0 * beam_energy * scattered_energy * sin2
    if q2 <= 0.0:
        return None
    w2 = mp2 + 2.0 * mp * energy_transfer - q2

    electron_energy_elastic = beam_energy * beam_energy * sin2
    electron_energy_elastic = electron_energy_elastic / (8.0 / 1.00797 * mp + 2.0 * beam_energy)
    scattered_energy_nuclear_elastic = beam_energy - electron_energy_elastic
    # This quantity is not used directly in the cross-section formula, but the
    # Fortran writes it out as (nu - nuel), which equals the excitation energy.
    excitation_energy = energy_transfer - electron_energy_elastic

    effective_potential = vcoul(a, z) if coulomb_correction else 0.0
    focusing_factor = 1.0 + effective_potential / beam_energy if coulomb_correction else 1.0
    beam_energy_effective = beam_energy + effective_potential
    scattered_energy_effective = scattered_energy + effective_potential

    q2_effective = 4.0 * beam_energy_effective * scattered_energy_effective * sin2
    if q2_effective <= 0.0:
        return None
    epsilon_effective = 1.0 / (1.0 + 2.0 * (energy_transfer * energy_transfer + q2_effective) / q2_effective * tan2)
    w2_effective = mp2 + 2.0 * mp * energy_transfer - q2_effective
    if abs(w2_effective - mp2) <= 0.0:
        return None
    bjorken_x_effective = q2_effective / (w2_effective - mp2 + q2_effective)
    kappa_effective = abs(w2_effective - mp2) / 2.0 / mp
    flux_effective = alpha * kappa_effective / (2.0 * pi2 * q2_effective) * scattered_energy_effective / beam_energy_effective / (1.0 - epsilon_effective)

    contribution_map: dict[str, float] = {}
    for kind, label in [(1, 'total_non_narrow'), (2, 'quasi_elastic'), (3, 'inelastic'), (4, 'mec')]:
        sigt, sigl = csfitcomp(w2_effective, q2_effective, a, z, xvalc, kind)
        sigma = flux_effective * (sigt + epsilon_effective * sigl)
        sigma = 0.3894e3 * 8.0 * pi2 * alpha / abs(w2_effective - mp2) * sigma
        sigma = sigma * focusing_factor * focusing_factor
        sigma = sigma / a
        contribution_map[label] = sigma

    narrow_state_total = 0.0
    for state in range(2, 22):
        narrow_state_total += nuccs12cs(z, a, beam_energy_effective, scattered_energy_effective, scattering_angle_deg, state) / 1000.0
    narrow_state_total = narrow_state_total / a

    cross_section_total = contribution_map['total_non_narrow'] + narrow_state_total

    # Follow the original Fortran printing cut.
    if not (scattered_energy > 0.01 and w2 < 40.0):
        return None

    return {
        'e': beam_energy,
        'theta': scattering_angle_deg,
        'nu': energy_transfer,
        'ep': scattered_energy,
        'ex': excitation_energy,
        'q2': q2,
        'w2': w2,
        'q2_coulomb': q2_effective,
        'w2_coulomb': w2_effective,
        'epsilon_coulomb': epsilon_effective,
        'flux_coulomb': flux_effective,
        'x_coulomb': bjorken_x_effective,
        'effective_potential': effective_potential,
        'focusing_factor': focusing_factor,
        'xs_total': cross_section_total,
        'xs_qe': contribution_map['quasi_elastic'],
        'xs_inelastic': contribution_map['inelastic'],
        'xs_mec': contribution_map['mec'],
        'xs_narrow_states': narrow_state_total,
        'xs_non_nuclear': contribution_map['total_non_narrow'],
    }


CROSS_SECTION_COLUMNS = [
    'e', 'theta', 'nu', 'ep', 'ex', 'q2', 'w2',
    'q2_coulomb', 'w2_coulomb', 'epsilon_coulomb', 'flux_coulomb',
    'x_coulomb', 'effective_potential', 'focusing_factor',
    'xs_total', 'xs_qe', 'xs_inelastic', 'xs_mec',
    'xs_narrow_states', 'xs_non_nuclear',
]



def calculate_cross_section_table(
    table: pd.DataFrame | np.ndarray | Iterable[tuple[float, float, float]],
    *,
    a: float = 12.0,
    z: float = 6.0,
    xvalc: Sequence[float] = XVALC,
    coulomb_correction: bool = True,
) -> pd.DataFrame:
    """Evaluate the qemodplot_norm-style cross section on a table of points.

    Parameters
    ----------
    table : DataFrame, ndarray, or iterable of tuples
        Input kinematics. Accepted forms are:

        - pandas DataFrame with columns ``e``, ``theta``, ``nu``
        - numpy array with shape ``(N, 3)`` storing ``[e, theta, nu]``
        - iterable of ``(e, theta, nu)`` tuples
    a : float, optional
        Nuclear mass number A. Default is 12.
    z : float, optional
        Nuclear charge Z. Default is 6.
    xvalc : sequence of float, optional
        1-based fit parameter table used by `csfitcomp`.
    coulomb_correction : bool, optional
        If True, apply the same Coulomb shift used by the Fortran driver.

    Returns
    -------
    pandas.DataFrame
        One row per valid kinematic point, including total cross section and
        component pieces.
    """
    if isinstance(table, pd.DataFrame):
        if not {'e', 'theta', 'nu'}.issubset(table.columns):
            raise ValueError("DataFrame input must contain columns: 'e', 'theta', 'nu'.")
        kinematics = table[['e', 'theta', 'nu']].to_numpy(dtype=float)
    elif isinstance(table, np.ndarray):
        arr = np.asarray(table, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError('NumPy input must have shape (N, 3) with columns [e, theta, nu].')
        kinematics = arr
    else:
        kinematics = np.asarray(list(table), dtype=float)
        if kinematics.ndim != 2 or kinematics.shape[1] != 3:
            raise ValueError('Iterable input must contain (e, theta, nu) triples.')

    rows = []
    for beam_energy, scattering_angle_deg, energy_transfer in kinematics:
        row = calculate_cross_section_point(
            float(beam_energy),
            float(scattering_angle_deg),
            float(energy_transfer),
            a=a,
            z=z,
            xvalc=xvalc,
            coulomb_correction=coulomb_correction,
        )
        if row is not None:
            rows.append(row)

    return pd.DataFrame(rows, columns=CROSS_SECTION_COLUMNS)
