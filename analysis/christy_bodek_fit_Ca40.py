"""Python translation of the Christy-Bodek universal fit Fortran code that
calculates fit values for RL RT response functions and cross-section of
electron-nucleus scattering.

Fortran core source codes:
    responseq.f
    qemodplot.f
    csfitcomp.f
    gsmearing.f
    qenuc21off.f
    mec2021.f
    sf.f
    rescsp.f / rescsn.f / resmodp.f / resmodn.f
    nuc12sf.f / nucffs12c.f / nucffs12ct.f
    formfacts.f

response output columns:
    qv, q2, ex, nu, rttot, rltot, rtqe, rlqe, rtie, rlie, rte, rle, rtns, rlns

cross-section output columns:
    e, theta, nu, ep, ex, q2, w2, q2_coulomb, w2_coulomb, epsilon_coulomb,
    flux_coulomb, x_coulomb, effective_potential, focusing_factor, xs_total,
    xs_qe, xs_inelastic, xs_mec, xs_narrow_state, xs_non_nuclear


Notes
-----
- use 1-based arrays for the fit parameter vectors (`xvalc`, resonance
  parameters, etc.) so the Python indexing matches the Fortran indexing.
- The algebra is kept close to the original code.
- A few numerical guards are added to avoid Python division-by-zero or
  sqrt-of-negative crashes at pathological kinematic points.
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
# Constants from response_Qvedges.f / other source files
# -----------------------------------------------------------------------------

MP_MAIN = 0.938273
MP = 0.938272
MN = 0.939565
ALPHA = 1.0 / 137.036
PI = 3.141593
PI2 = PI * PI


XVALC = _pad1([
    0.91648e-01, 0.12714e+02, 0.13380e+00, 0.69068e+01, 0.77023e+00,
    0.76437e-01, 0.87115e+01, 0.18976e+01, 0.66472e+00, -0.39215e+01,
    0.99320e+00, 0.98312e+00, 0.10302e+01, 0.10009e+01, 0.10000e+01,
    0.10070e+01, 0.97472e+00, 0.10059e+01, 0.98892e+00, 0.99434e+00,
    0.10000e+01, 0.99596e+00, 0.10028e+01, 0.10122e+01, 0.10045e+01,
    0.79845e+00, 0.11295e-05, -0.97071e+00, 0.92502e+00, 0.20146e+01,
    # 0.24416e+01, 0.24499e+01, 0.31154e+01, 0.72998e+00, 0.22800e+00,
    # 0.76502e-02, 0.25718e+00, 0.31429e-01, 0.58780e-01, -0.15059e+00,
    # FIXME: xvlac[35], xvlac[37] for Ca40
    0.24416e+01, 0.24499e+01, 0.31154e+01, 0.72998e+00, 0.26000e+00,
    0.76502e-02, 0.29000e+00, 0.31429e-01, 0.58780e-01, -0.15059e+00,
    0.38790e-01, 0.77051e-01, 0.26795e+00, 0.17673e+00, 0.10451e-01,
])


RESCSP_XVAL = _pad1([
    0.12291e+01, 0.15173e+01, 0.15044e+01, 0.17100e+01, 0.16801e+01,
    0.14312e+01, 0.12616e+00, 0.23000e+00, 0.92594e-01, 0.90606e-01,
    0.75000e-01, 0.35067e+00, 0.75729e+01, 0.56091e+01, 0.94606e+01,
    0.20156e+01, 0.66190e+01, 0.41732e+00, 0.23980e-01, 0.53136e+01,
    0.63752e+00, 0.11484e+02, 0.69949e-01, 0.26191e+01, 0.53603e-01,
    0.65000e+02, 0.15351e+00, 0.20624e+01, 0.23408e+01, 0.16100e+02,
    0.62414e+02, 0.17201e+01, 0.23261e+00, 0.65000e+02, 0.23292e+01,
    0.14980e+01, 0.23000e+00, 0.63385e+00, 0.19093e-01, 0.61061e-01,
    0.29146e-02, 0.54388e+00, 0.77997e+00, 0.28783e+00, 0.10605e+01,
    0.69793e+00, 0.20009e+01, 0.57000e+00, 0.41632e+01, 0.38427e+00,
    0.10000e+01, 0.99842e+00, 0.98719e+00, 0.10168e+01, 0.98945e+00,
    0.99594e+00, 0.98799e+00, 0.10271e+01, 0.10650e+01, 0.97920e+00,
    0.10152e+01, 0.99622e+00, 0.81011e+01, 0.10070e-02, 0.14857e+01,
    0.33445e+01, 0.31641e-09, 0.69755e+02, 0.55228e+01, 0.14438e+00,
    0.60474e+01, 0.65395e-07, 0.14129e+01, 0.58609e+00, 0.36220e+01,
    0.92699e+00, 0.14418e+01, 0.86403e-02, 0.10001e-03, 0.75106e+00,
    0.76077e+00, 0.42272e+00, 0.55511e-11, 0.52486e+00, 0.58153e+00,
    0.15798e+01, 0.50105e+00, 0.89149e+02, 0.72789e+00, 0.24813e-01,
    -0.61906e+00, 0.10000e+01, 0.00000e+00, 0.00000e+00, 0.68158e+03,
    0.12429e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.10000e-05,
])


RESCSN_XVAL = _pad1([
    0.12291e+01, 0.15173e+01, 0.15044e+01, 0.17100e+01, 0.16801e+01,
    0.14312e+01, 0.12616e+00, 0.23000e+00, 0.92594e-01, 0.90606e-01,
    0.75000e-01, 0.35067e+00, 0.69500e+01, 0.86607e+01, 0.11555e+02,
    0.22138e+01, 0.44887e+01, 0.11332e+03, 0.46513e+03, 0.31173e+01,
    0.96301e+00, 0.14956e+00, 0.20761e-07, 0.10440e+01, 0.40143e-03,
    0.90028e+02, 0.75248e-01, 0.20532e+00, 0.14359e-01, 0.29830e+03,
    0.19935e+00, 0.26923e+01, 0.48616e+01, 0.86000e+02, 0.67813e+04,
    0.44281e+02, 0.29572e+00, 0.65424e+00, 0.23787e-09, 0.52052e-01,
    0.39926e-08, 0.29941e+00, 0.97516e+00, 0.46934e-01, 0.14246e+03,
    0.55801e+00, 0.19349e+01, 0.27400e+00, 0.38891e+00, 0.40000e-02,
    0.10108e+01, 0.97020e+00, 0.98248e+00, 0.97768e+00, 0.10425e+01,
    0.10198e+01, 0.97822e+00, 0.98239e+00, 0.10103e+01, 0.10076e+01,
    0.10044e+01, 0.99687e+00, 0.16696e+01, 0.10721e-06, 0.54114e+00,
    0.79999e+03, 0.56715e+02, 0.73797e+03, 0.37949e+02, 0.14796e+03,
    0.30498e+01, 0.24459e+00, 0.95574e+00, 0.35577e+00, 0.21228e-05,
    0.96696e+01, 0.27563e+01, 0.93024e-01, 0.33559e+02, 0.31207e-01,
    0.29020e+02, 0.86417e+00, 0.36471e-08, 0.99167e+00, 0.68124e+00,
    0.10000e-01, 0.90227e-01, 0.40115e+01, 0.29915e+01, 0.45929e-01,
    -0.16758e+01, 0.78493e+01, 0.78184e+01, 0.42074e+01, 0.41179e-05,
    0.80597e+00, 0.00000e+00, 0.00000e+00, 0.10045e+01, 0.62364e+00,
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
# resmodp.f / resmodn.f helpers
# -----------------------------------------------------------------------------


def _resmod_common(sf: int, w2: float, q2: float, xval: Sequence[float], *, proton: bool) -> float:
    """Shared implementation of RESMODP / RESMODN.

    The algebra follows the Fortran very closely, with only tiny numerical
    guards to keep Python well-defined at thresholds.
    """
    mp = 0.9382727 if proton else 0.939565
    mpi = 0.134977
    meta = 0.547862
    mp2 = mp * mp
    w = math.sqrt(max(w2, 0.0))
    wdif1 = w - (mp + mpi)
    wdif2 = w - (mp + meta)
    q20 = xval[50]

    # Branching ratios and orbital angular momenta.
    br = [[0.0] * 4 for _ in range(8)]  # use [1..7][1..3]
    ang = [0.0] * 8
    for i, v in enumerate([1.0, 0.45, 0.60, 0.65, 0.60, 0.65, 0.60], start=1):
        br[i][1] = v
    for i, v in enumerate([0.0, 0.40, 0.08, 0.0, 0.20, 0.0, 0.0], start=1):
        br[i][3] = v
    for i in range(1, 8):
        br[i][2] = 1.0 - br[i][1] - br[i][3]
    for i, v in enumerate([1.0, 0.0, 2.0, 3.0, 0.0, 1.0, 3.0], start=1):
        ang[i] = v

    x0 = [0.0] * 8
    for i in range(1, 8):
        x0[i] = 0.160 if proton else 0.16
    if sf == 2:
        x0[1] = 0.07

    mon = 1.0 / (1.0 + q2 / 1.5)

    xb = q2 / (q2 + w2 - mp2)
    xpr1 = 1.0 + (w2 - (mp + mpi) ** 2) / (q2 + q20)
    xpr1 = 1.0 / xpr1
    xpr2 = 1.0 + (w2 - (mp + meta) ** 2) / (q2 + q20)
    xpr2 = 1.0 / xpr2
    if w <= (mp + mpi):
        xpr1 = 1.0
    if w <= (mp + meta):
        xpr2 = 1.0

    # Threshold kinematics for Breit-Wigner factors.
    k = (w2 - mp2) / (2.0 * mp)
    kcm = (w2 - mp2) / (2.0 * w) if w > 0.0 else 0.0

    epicm = (w2 + mpi**2 - mp2) / (2.0 * w) if w > 0.0 else 0.0
    ppicm = _sqrt_pos(epicm**2 - mpi**2)
    epi2cm = (w2 + (2.0 * mpi) ** 2 - mp2) / (2.0 * w) if w > 0.0 else 0.0
    ppi2cm = _sqrt_pos(epi2cm**2 - (2.0 * mpi) ** 2)
    eetacm = (w2 + meta * meta - mp2) / (2.0 * w) if w > 0.0 else 0.0
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
    ppicmr = [0.0] * 8
    ppi2cmr = [0.0] * 8
    petacmr = [0.0] * 8
    pgam = [0.0] * 8
    pwid = [[0.0] * 4 for _ in range(8)]

    for i in range(1, 8):
        mi = mass[i]
        mi2 = mi * mi
        kr[i] = (mi2 - mp2) / (2.0 * mp)
        kcmr[i] = (mi2 - mp2) / (2.0 * mi)
        epicmr = (mi2 + mpi**2 - mp2) / (2.0 * mi)
        ppicmr[i] = _sqrt_pos(epicmr**2 - mpi**2)
        epi2cmr = (mi2 + (2.0 * mpi) ** 2 - mp2) / (2.0 * mi)
        ppi2cmr[i] = _sqrt_pos(epi2cmr**2 - (2.0 * mpi) ** 2)
        eetacmr = (mi2 + meta * meta - mp2) / (2.0 * mi)
        petacmr[i] = _sqrt_pos(eetacmr**2 - meta**2)

        # Partial widths.
        if ppicmr[i] > 0.0:
            pwid[i][1] = intwidth[i] * (ppicm / ppicmr[i]) ** (2.0 * ang[i] + 1.0) * (
                (ppicmr[i] ** 2 + x0[i] ** 2) / (ppicm**2 + x0[i] ** 2)
            ) ** ang[i]
        else:
            pwid[i][1] = 0.0

        if ppi2cmr[i] > 0.0:
            pwid[i][2] = intwidth[i] * (ppi2cm / ppi2cmr[i]) ** (2.0 * ang[i] + 4.0) * (
                (ppi2cmr[i] ** 2 + x0[i] ** 2) / (ppi2cm**2 + x0[i] ** 2)
            ) ** (ang[i] + 2.0)
            pwid[i][2] = (w / mi) * pwid[i][2]
        else:
            pwid[i][2] = 0.0

        pwid[i][3] = 0.0
        if i in (2, 5) and petacmr[i] > 0.0:
            pwid[i][3] = intwidth[i] * (petacm / petacmr[i]) ** (2.0 * ang[i] + 1.0) * (
                (petacmr[i] ** 2 + x0[i] ** 2) / (petacm**2 + x0[i] ** 2)
            ) ** ang[i]

        if kcmr[i] > 0.0:
            pgam_factor = (kcm / kcmr[i]) ** 2 * (kcmr[i] ** 2 + x0[i] ** 2) / (kcm**2 + x0[i] ** 2)
            pgam[i] = intwidth[i] * pgam_factor
        else:
            pgam[i] = 0.0

        width[i] = br[i][1] * pwid[i][1] + br[i][2] * pwid[i][2] + br[i][3] * pwid[i][3]

    # Q^2 dependence of resonance heights.
    height = [0.0] * 8
    rescoef = [[0.0] * 5 for _ in range(7)]
    for i in range(1, 7):
        for j in range(1, 5):
            num += 1
            rescoef[i][j] = xval[num]
        if sf == 1:
            height[i] = rescoef[i][1] * (
                1.0 + rescoef[i][2] * q2 / (1.0 + rescoef[i][3] * q2)
            ) * mon ** rescoef[i][4]
        else:
            height[i] = (rescoef[i][1] + rescoef[i][2] * q2) * _exp(-1.0 * rescoef[i][3] * q2)
        height[i] = height[i] * height[i]

    if proton:
        if sf == 2:
            height[7] = (xval[16] + xval[20] * q2) * _exp(-1.0 * xval[24] * q2)
        else:
            height[7] = xval[49] * mon ** xval[45]
    else:
        if sf == 2:
            height[7] = (xval[44] + xval[45] * q2) * _exp(-1.0 * xval[46] * q2)
        else:
            height[7] = xval[49] * mon
    height[7] = height[7] * height[7]

    nr_coef = [[0.0] * 5 for _ in range(4)]
    for i in range(1, 4):
        for j in range(1, 5):
            num += 1
            nr_coef[i][j] = xval[num]

    # Breit-Wigner sum.
    sig_res = 0.0
    if abs(k) > 0.0 and abs(kcm) > 0.0:
        for i in range(1, 8):
            denom = (w2 - mass[i] ** 2) ** 2 + (mass[i] * width[i]) ** 2
            if denom <= 0.0 or intwidth[i] == 0.0:
                sigr = 0.0
            else:
                sigr = width[i] * pgam[i] / denom
                sigr = height[i] * kr[i] / k * kcmr[i] / kcm * sigr / intwidth[i]
            sig_res += sigr
    sig_res *= w
    if sf == 2:
        sig_res *= q2

    # Non-resonant background.
    sig_nr = 0.0
    if sf == 1 and xpr1 < 1.0:
        if proton:
            a0 = xval[37] / (1.0 + q2 / xval[42]) ** xval[43]
            t1 = xval[38] * math.log(1.06 + q2) + xval[39] / math.log(1.06 + q2)
            t2 = xval[40] * (1.0 + q2 / xval[41]) ** xval[44]
        else:
            a0 = xval[37] / (1.0 + q2 / xval[42]) ** xval[43]
            t1 = xval[38] * math.log(1.05 + q2) + xval[39] / (1.05 + q2)
            t2 = xval[40] * (1.0 + q2 / xval[41]) ** xval[44]

        if xpr1 <= 1.0:
            sig_nr = 389.4 * a0 * (1.0 - xpr1) ** t1 * xpr1**t2
        if xpr2 <= 1.0:
            sig_nr += xval[46] * 389.4 * a0 * (1.0 - xpr2) ** t1 * xpr2**t2

    elif sf == 2 and xpr1 < 1.0:
        if proton:
            a0 = xval[37] / (1.0 + q2 / xval[39]) ** 2.0
            t1 = xval[38] / (1.0 + q2 / xval[40]) + xval[32] * math.log(q2 + xval[36])
            t2 = xval[41]
        else:
            a0 = xval[37] / (1.0 + q2 / xval[39]) ** 2.0
            t1 = xval[38] / (1.0 + q2 / xval[40]) + xval[32] * math.log(q2 + xval[36])
            t2 = xval[41] / (1.0 + q2 / xval[42]) ** xval[43]

        if xpr1 <= 1.0:
            sig_nr += 389.4 * a0 * xb * (1.0 - xpr1) ** t1 * xpr1**t2

    sig = sig_res + sig_nr
    if (w - mp) < wdif1:
        sig = 0.0
    _ = wdif2  # kept for parity with Fortran; not used further.
    _ = nr_coef  # read for structural closeness.
    return sig


# -----------------------------------------------------------------------------
# rescsp.f / rescsn.f / sf.f
# -----------------------------------------------------------------------------


def rescsp(w2: float, q2: float) -> Tuple[float, float]:
    """rescsp(W2,Q2,sigT,sigL)."""
    xval1 = [0.0] * 51
    xvalL = [0.0] * 51
    for i in range(1, 51):
        xval1[i] = RESCSP_XVAL[i]
        xvalL[i] = RESCSP_XVAL[50 + i]
        if i <= 12:
            xvalL[i] = xval1[i]
        if i in (47, 48):
            xvalL[i] = xval1[i]
    sigt = _resmod_common(1, w2, q2, xval1, proton=True)
    sigl = _resmod_common(2, w2, q2, xvalL, proton=True)
    return sigt, sigl


def rescsn(w2: float, q2: float) -> Tuple[float, float]:
    """rescsn(W2,Q2,sigtn,sigln)."""
    xval1 = [0.0] * 51
    xvalL = [0.0] * 51
    for i in range(1, 51):
        xval1[i] = RESCSN_XVAL[i]
        xvalL[i] = RESCSN_XVAL[50 + i]
        if i <= 12:
            xvalL[i] = xval1[i]
        if i in (47, 48):
            xvalL[i] = xval1[i]
    sigtn = _resmod_common(1, w2, q2, xval1, proton=False)
    sigln = _resmod_common(2, w2, q2, xvalL, proton=False)
    return sigtn, sigln


def sf(w2: float, q2: float) -> Tuple[float, float, float, float, float, float]:
    """SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)."""
    mp = 0.938272
    mp2 = mp * mp
    alpha = 1.0 / 137.03599
    x = q2 / (q2 + w2 - mp2)

    sigtp, siglp = rescsp(w2, q2)
    sigtn, sigln = rescsn(w2, q2)

    pref = abs(w2 - mp2) / (0.3894e3 * PI2 * alpha * 8.0)
    f1p = sigtp * pref
    f1n = sigtn * pref
    flp = siglp * 2.0 * x * pref
    fln = sigln * 2.0 * x * pref
    denom = 1.0 + 4.0 * mp2 * x * x / q2
    f2p = (2.0 * x * f1p + flp) / denom
    f2n = (2.0 * x * f1n + fln) / denom
    return f1p, flp, f2p, f1n, fln, f2n


# -----------------------------------------------------------------------------
# gsmearing.f
# -----------------------------------------------------------------------------


def gsmearing(z: float, a: float, w2: float, q2: float, xvalc: Sequence[float]) -> Tuple[float, float, float]:
    """GSMEARING(Z, A, W2, Q2, xvalc, F1, F2, FL)."""
    nbins = 98
    nwid = 3.3
    bw = 2.0 * nwid / float(nbins)

    exmin = 0.0165
    mp = 0.938272
    mp2 = mp * mp
    x = q2 / (q2 + w2 - mp2)
    nu = (w2 + q2 - mp2) / (2.0 * mp)
    qv = math.sqrt(nu * nu + q2)
    kappa2 = 1.0 + 4.0 * mp2 * x * x / q2
    nuel = q2 / (2.0 * (0.931494 * a))
    ex = nu - nuel

    es = 0.008
    if a >= 3.0:
        kf = xvalc[37]
        qvt = min(qv, 1.0)
        es = xvalc[38]
        es = es - xvalc[39] * (1.0 - qvt)

    norm = math.sqrt(PI)
    ncor = 1.000
    norm = norm / ncor

    f1p = f1n = f2p = f2n = flp = fln = 0.0
    fytot = 0.0
    fytot2 = 0.0

    pf = 0.5 * kf
    pf2 = pf * 1.5
    dw2dpf = 2.0 * qv
    dw2des = 2.0 * (nu + mp)

    for ism in range(1, nbins + 1):
        xxp = -nwid + bw * float(ism - 1)
        fyuse = bw / math.sqrt(2.0) / norm * _exp(-0.5 * xxp * xxp)

        wsqp = w2 + xxp * pf * dw2dpf - es * dw2des
        wsqp2 = w2 + xxp * pf2 * dw2dpf - es * dw2des

        fytot += fyuse
        fytot2 += fyuse
        _ = fytot, fytot2, wsqp2  # kept for structural parity

        frac = 0.0
        for j in range(1, 2):
            if j == 1:
                fract = 1.0 - frac
                w2t = wsqp
            else:
                fract = frac
                w2t = wsqp2

            if w2t > 1.159:
                xt = q2 / (q2 + w2t - mp2)
                xp = 1.0 + (w2t - mp) / (q2 + xvalc[34])
                xp = 1.0 / xp
                offshell = 1.0

                emcfac = (xvalc[26] + xvalc[27] * xp * xp) / (1.0 + xvalc[28] * xp + xvalc[29] * xp * xp)
                emcfacL = xvalc[30] * (1.0 + xvalc[31] * xp * xp) * (1.0 + xvalc[32] * xp * xp) * _exp(-1.0 * xvalc[33] * xp)

                f1pp, flpp, f2pp, f1nn, flnn, f2nn = sf(w2t, q2)

                f1pp *= emcfac * offshell
                f1nn *= emcfac * offshell
                flpp *= emcfac * emcfacL * offshell
                flnn *= emcfac * emcfacL * offshell
                denom = 1.0 + 4.0 * xt * xt * mp2 / q2
                f2pp = (2.0 * xt * f1pp + flpp) / denom
                f2nn = (2.0 * xt * f1nn + flnn) / denom

                f1p += f1pp * fyuse * fract
                f1n += f1nn * fyuse * fract
                f2p += f2pp * fyuse * fract
                f2n += f2nn * fyuse * fract
                flp += flpp * fyuse * fract
                fln += flnn * fyuse * fract

    f2 = z * f2p + (a - z) * f2n
    fl = z * flp + (a - z) * fln

    cof = math.sqrt(max(ex - exmin, 0.0)) / math.sqrt(0.025 - exmin)
    cof = _clamp(cof, 0.0, 1.0)
    if ex <= exmin:
        cof = 0.0
    f2 *= cof
    fl *= cof

    f1 = (kappa2 * f2 - fl) / (2.0 * x)
    f1 = max(f1, 0.0)
    f2 = max(f2, 0.0)
    fl = max(fl, 0.0)
    return f1, f2, fl


# -----------------------------------------------------------------------------
# qenuc21off.f
# -----------------------------------------------------------------------------


def qenuc21off(z: float, a: float, q2: float, w2: float, xvalc: Sequence[float]) -> Tuple[float, float]:
    """QENUC21OFF(Z, A, Q2, W2, xvalc, F1, F2)."""
    mp = 0.938272
    psimax = 5.0
    moff = 1.0 * mp
    paulitype = 1

    f1 = 0.0
    f2 = 0.0
    ia = int(a)
    avgn = a - z
    if ia == 1:
        return 0.0, 0.0

    nuel = q2 / (2.0 * (0.931494 * a))
    nu = (w2 - mp**2 + q2) / (2.0 * mp)
    if nu <= 0.0 or q2 <= 0.0:
        return 0.0, 0.0

    tau = q2 / (4.0 * mp**2)
    qv = math.sqrt(nu * nu + q2)
    ex = nu - nuel

    gmp, gep, gmn, gen = formfacts(q2)

    if ia > 2:
        esmin = 0.020
        if ia == 12:
            # esmin = 0.0165
            esmin = 0.022 # temporary
        kf = xvalc[35]
        es = esmin + xvalc[36]
        if qv < 1.5:
            es = es - xvalc[36] * (1.0 - qv / 1.5) ** 0.1
    else:
        # The original code only targets A>2 here; keep a safe fallback.
        esmin = 0.020
        kf = xvalc[35]
        es = esmin

    nup = nu - es

    if (qv > 2.0 * kf) or (ia == 1):
        pauli_sup2 = 1.0
    else:
        pauli_sup2 = 0.75 * (qv / kf) * (1.0 - ((qv / kf) ** 2) / 12.0)
    pauli_sup1 = pauli_sup2
    _ = pauli_sup1

    kappa = qv / (2.0 * mp)
    lam = nu / (2.0 * mp)
    lamp = nup / (2.0 * mp)
    lampn = -lamp
    taup = kappa**2 - lamp**2
    xi = math.sqrt(1.0 + (kf / moff) ** 2) - 1.0

    # The Fortran assumes the kinematics are in a valid QE region.
    if taup <= 0.0 or (1.0 + taup) <= 0.0:
        return 0.0, 0.0

    root_common = math.sqrt(taup * (1.0 + taup))
    den_psi = math.sqrt(xi) * math.sqrt((1.0 + lam) * taup + kappa * root_common)
    den_psip = math.sqrt(xi) * math.sqrt((1.0 + lamp) * taup + kappa * root_common)
    den_psipn_arg = (1.0 + lampn) * taup + kappa * root_common
    if den_psi == 0.0 or den_psip == 0.0 or den_psipn_arg <= 0.0:
        return 0.0, 0.0

    psi = (lam - taup) / den_psi
    psip = (lamp - taup) / den_psip
    psipn = (lampn - taup) / (math.sqrt(xi) * math.sqrt(den_psipn_arg))

    nuL = (q2 / qv / qv) ** 2
    nuT = tau / (2.0 * kappa**2)

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

    num = 2.0 * kappa * (1.0 + xi * (1.0 + psi**2) / 2.0)
    gl /= num
    gt /= num

    fy = 1.5576 / (1.0 + 1.7720**2 * (psip + 0.3014) ** 2) / (1.0 + _exp(-2.4291 * psip))
    fyp = xvalc[7] / (1.0 + xvalc[8] ** 2 * (psip + xvalc[9]) ** 2) / (1.0 + _exp(xvalc[10] * psip)) * (1.0 - abs(psip) / psimax) ** 2.0 * (1.0 + abs(psip) / psimax) ** 2.0
    fyn = xvalc[7] / (1.0 + xvalc[8] ** 2 * (psipn + xvalc[9]) ** 2) / (1.0 + _exp(xvalc[10] * psipn)) * (1.0 - abs(psipn) / psimax) ** 2.0 * (1.0 + abs(psipn) / psimax) ** 2.0

    if psip > psimax:
        fyp = 0.0
    if psipn > psimax:
        fyn = 0.0
    fyp = max(0.0, fyp)
    fyn = max(0.0, fyn)

    if paulitype == 1:
        fy = max(0.0, fyp - fyn)
        x2 = qv / kf
        pb2L = 1.0 - xvalc[40] * (4.0 - x2) ** 2.5 - xvalc[41] * (4.0 - x2) ** 3.5
        pb2L = pb2L - xvalc[44] * (4.0 - x2) ** 1.5
        pb2L = pb2L * (x2 - 0.2) ** 2 / (x2 - 0.18) ** 2.0
        if x2 > 4.0:
            pb2L = 1.0
        pb2L = _clamp(pb2L, 0.0, 1.0)
        if x2 > 4.0:
            pb2L = 1.0
        pb2L = min(pb2L, 1.0)
        if psip > psimax:
            fy = 0.0
        fyl = pb2L * fy
    else:
        fy = pauli_sup2 * fyp
        fyl = fy

    f2 = nu / kf * (fyl * nuL * gl + fy * nuT * gt)
    f1 = mp * fy / kf * gt / 2.0

    cof = math.sqrt(max(ex - esmin, 0.0)) / math.sqrt(0.025 - esmin)
    # cof = math.sqrt(max(ex - esmin, 0.0)) / math.sqrt(0.038 - esmin) # FIXME: temporary

    cof = _clamp(cof, 0.0, 1.0)
    if ex <= esmin:
        cof = 0.0
    f1 *= cof
    f2 *= cof
    f2 = max(f2, 0.0)
    f1 = max(f1, 0.0)
    return f1, f2


# -----------------------------------------------------------------------------
# mec2021.f
# -----------------------------------------------------------------------------


def mec2021(z: float, a: float, w2: float, q2: float, xvalm: Sequence[float]) -> float:
    """MEC2021(z,a,w2,q2,xvalm,f1mec)."""
    mp = 0.938272
    q20 = 0.00001
    mp2 = mp * mp

    f1mec = 0.0
    if w2 <= 0.0:
        return 0.0

    w = math.sqrt(w2)
    nu = (w2 - mp2 + q2) / (2.0 * mp)
    x = q2 / (2.0 * mp * nu)
    qv2 = q2 + nu**2
    _ = w, qv2
    numin = 0.0165
    w2min = mp2 + 2.0 * mp * numin - q2
    xmax = q2 / (2.0 * mp * numin)

    nuel = q2 / (2.0 * (0.931494 * a))
    ex = nu - nuel

    if a < 2.5:
        return 0.0

    a1 = xvalm[1]
    a2 = xvalm[45]

    y = a * _exp(-1.0 * q2 * q2 / xvalm[2]) * (q2 + q20) ** 2 / (xvalm[3] + q2) ** xvalm[4]
    b1 = xvalm[5]
    b2 = 1.275
    c1 = xvalm[42] + xvalm[43] * q2
    c2 = 0.265

    t1 = (w2 - b1) ** 2 / (2.0 * c1**2)
    t2 = (w2 - b2) ** 2 / (2.0 * c2**2)
    dw2 = w2 - w2min
    if dw2 < 0.0:
        dw2 = 0.0

    f1mec = y * _exp(-1.0 * t1)
    f1mec = f1mec * (dw2) ** 1.5
    f1mec2 = a2 * y * (q2 + q20) ** 1.5 * _exp(-1.0 * t2)
    f1mec = a1 * f1mec + f1mec2

    if nu < numin:
        f1mec = 0.0
    cof = math.sqrt(max(ex - numin, 0.0)) / math.sqrt(0.025 - numin)
    cof = _clamp(cof, 0.0, 1.0)
    if ex <= numin:
        cof = 0.0
    f1mec *= cof

    if dw2 <= 0.0 or x > xmax:
        f1mec = 0.0
    if f1mec <= 1.0e-9 or x >= xmax:
        f1mec = 0.0
    return f1mec


# -----------------------------------------------------------------------------
# nucffs12c.f / nucffs12ct.f / nuc12sf.f
# -----------------------------------------------------------------------------


def nucffs12c(a: float, z: float, q32: float, state: int) -> float:
    """NUCFFS12C(A,Z,q32,state,FF). Only the 12C states used by nuc12sf are needed."""
    q2f = q32 / (0.1975 * 0.1975)
    qf = math.sqrt(max(q2f, 0.0))
    radius = 2.45
    ff2 = 0.0

    if state == 1:
        # In the original Fortran, Q2 is undeclared/unset here. Using q32 is the
        # physically sensible interpretation, but this branch is not used by the
        # current driver (which sums states 2..22).
        x2 = (q32 / 0.197328**2) * radius**2
        alp = (z - 2.0) / 3.0
        char = x2 * (2.0 + 3.0 * alp) / (12.0 + 30.0 * alp)
        ff = 0.0
        if char < 80.0:
            ff = _exp(-char) * (1.0 - alp * x2 / (6.0 + 15.0 * alp))
        ff2 = ff * ff
        alph = 4.0 / 3.0
        a0 = 1.65
        g0 = 8.0e-5 * _exp(-1.0 * ((q2f - 2.9) / 0.44) ** 2.0)
        h0 = (1.0 - alph * q2f * a0 * a0 / 2.0 / (2.0 + 3.0 * alph)) * _exp(-q2f * a0 * a0 / 4.0)
        g4 = 1.0e-5 * _exp(-1.0 * ((q2f - 4.0) / 1.2) ** 2.0)
        ff2 = h0 * h0 + g0 + g4
    elif state == 2:
        g1 = 1.41e-2 * _exp(-1.0 * (q2f - 1.125) ** 2 / 1.71 / 1.71)
        g2 = 7.0e-4 * _exp(-1.0 * (q2f - 3.7) ** 2 / 1.6 / 1.6)
        g3 = 3.3e-6 * _exp(-1.0 * (q2f - 6.5) ** 2 / 7.0 / 7.0)
        g4 = 0.0
        ff2 = q2f**3 / (q2f**3 + 0.1) * (g1 + g2 + g3 + g4)
    elif state == 3:
        b = 1.3457
        a1 = 0.52 * (b * qf) ** 2
        a2 = -0.025 * (b * qf) ** 4
        a3 = -0.7e-2 * (b * qf) ** 6
        a4 = 0.5e-3 * (b * qf) ** 8
        a5 = -0.5e-4 * (b * qf) ** 10
        ff2 = (1.0 / z) * _exp(-0.5 * b * b * q2f) * (a1 + a2 + a3 + a4 + a5)
        ff2 = ff2 * ff2
    elif state == 4:
        g1 = 5.0e-3 * _exp(-1.0 * (q2f - 1.46) ** 2 / 1.6 / 1.6)
        g2 = 6.6e-4 * _exp(-1.0 * (q2f - 3.46) ** 2 / 2.0 / 2.0)
        g3 = 7.0e-6 * _exp(-1.0 * (q2f - 7.0) ** 2 / 2.8 / 2.8)
        ff2 = q2f**3 / (q2f**3 + 0.2) * (g1 + g2 + g3)
    elif state == 5:
        g1 = 5.0e-4 * _exp(-1.0 * (qf - 1.0) ** 2 / 0.3 / 0.3)
        g2 = 8.0e-4 * _exp(-1.0 * (qf - 1.4) ** 2 / 0.4 / 0.4)
        g3 = 0.0 * _exp(-1.0 * (qf - 7.0) ** 2 / 2.5 / 2.5)
        e1 = 0.0
        ff2 = g1 + g2 + g3 + e1
    elif state == 6:
        g1 = 4.0e-4 * _exp(-1.0 * (qf - 1.0) ** 2 / 0.35 / 0.35)
        g2 = 8.0e-4 * _exp(-1.0 * (qf - 1.75) ** 2 / 0.45 / 0.45)
        g3 = 4.0e-4 * _exp(-1.0 * (qf - 0.85) ** 2 / 0.65 / 0.65)
        e1 = 0.0
        g4 = 0.0
        ff2 = g1 + g2 + g3 + e1 + g4
    elif state == 7:
        g1 = 6.0e-4 * _exp(-1.0 * (qf - 0.85) ** 2 / 0.7 / 0.7)
        ff2 = g1
    elif state == 8:
        g1 = 12.0e-4 * _exp(-1.0 * (qf - 1.05) ** 2 / 0.6 / 0.6)
        ff2 = g1
    elif state == 9:
        g1 = 3.2e-4 * _exp(-1.0 * (qf - 1.3) ** 2 / 0.5 / 0.5)
        ff2 = g1
    elif state == 10:
        g1 = 1.6e-4 * _exp(-1.0 * (qf - 1.2) ** 2 / 0.42 / 0.42)
        g2 = 1.6e-5 * _exp(-1.0 * (qf - 1.8) ** 2 / 0.4 / 0.4)
        ff2 = g1 + g2
    elif state == 11:
        g1 = 2.8e-3 * _exp(-1.0 * (qf - 0.60) ** 2 / 0.15 / 0.15)
        g2 = 6.9e-3 * _exp(-1.0 * (qf - 0.84) ** 2 / 0.55 / 0.55)
        ff2 = g1 + g2
    elif state == 12:
        g1 = 0.0047 * _exp(-1.0 * (qf - 1.0) ** 2 / 0.48 / 0.48)
        ff2 = g1
    elif state == 13:
        g1 = 0.0026 * _exp(-1.0 * (qf - 1.49) ** 2 / 0.7 / 0.7)
        ff2 = g1

    if qf * qf > 12.0:
        ff2 = 0.0
    ff2 = max(ff2, 0.0)
    return math.sqrt(ff2)



def nucffs12ct(a: float, z: float, q2: float, state: int) -> float:
    """NUCFFS12CT(A,Z,Q2,state,FF)."""
    q2f = q2 / (0.1975 * 0.1975)
    qf = math.sqrt(max(q2f, 0.0))
    ff2 = 0.0

    if state == 14:
        g1 = 2.5e-4 * _exp(-1.0 * (qf - 0.63) ** 2 / 0.4 / 0.4)
        g2 = 2.8e-4 * _exp(-1.0 * (qf - 0.84) ** 2 / 0.2 / 0.2)
        g3 = 0.0
        e1 = -2.5e-5 * _exp(-1.0 * qf)
        ff2 = g1 + g2 + g3 + e1
    elif state == 15:
        g1 = 5.9e-4 * _exp(-1.0 * (qf - 1.2) ** 2 / 0.55 / 0.55)
        g2 = 2.4e-4 * _exp(-1.0 * (qf - 2.2) ** 2 / 0.6 / 0.6)
        ff2 = g1 + g2
    elif state == 16:
        g1 = 2.6e-4 * _exp(-1.0 * (qf - 1.6) ** 2 / 0.6 / 0.6)
        g2 = 5.0e-5 * _exp(-1.0 * (qf - 2.5) ** 2 / 0.35 / 0.35)
        ff2 = g1 + g2
    elif state == 17:
        g1 = 2.1e-4 * _exp(-1.0 * (qf - 0.8) ** 2 / 0.35 / 0.35)
        g2 = 1.6e-4 * _exp(-1.0 * (qf - 1.2) ** 2 / 0.425 / 0.425)
        ff2 = g1 + g2
    elif state == 18:
        g1 = 9.5e-4 * _exp(-1.0 * (qf - 1.27) ** 2 / 0.77 / 0.77)
        g2 = 3.5e-4 * _exp(-1.0 * (qf - 1.7) ** 2 / 0.6 / 0.6)
        g3 = 1.0e-4 * _exp(-1.0 * (qf - 2.2) ** 2 / 0.3 / 0.3)
        e1 = -3.6e-4 * _exp(-1.0 * qf)
        ff2 = g1 + g2 + g3 + e1
    elif state == 19:
        g1 = 2.1e-4 * _exp(-1.0 * (qf - 1.45) ** 2 / 0.5 / 0.5)
        g2 = 5.5e-5 * _exp(-1.0 * (qf - 2.1) ** 2 / 0.4 / 0.4)
        ff2 = g1 + g2
    elif state == 20:
        g1 = 1.8e-3 * _exp(-1.0 * (qf - 0.8) ** 2 / 0.36 / 0.36)
        g2 = 1.0e-4 * _exp(-1.0 * (qf - 1.5) ** 2 / 0.5 / 0.5)
        ff2 = g1 + g2
    elif state == 21:
        g1 = 9.0e-4 * _exp(-1.0 * (qf - 0.35) ** 2 / 0.3 / 0.3)
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        ff2 = g1 + g2 + g3 + e1

    if qf * qf > 12.0:
        ff2 = 0.0
    ff2 = max(ff2, 0.0)
    return math.sqrt(ff2)



def nuc12sf(z: float, a: float, nu: float, q2p: float, state: int) -> Tuple[float, float]:
    """NUC12SF(Z,A,nu,q2p,state,F1,FL)."""
    mp = 0.93827
    pi = 3.14159

    # NOTE: the original Fortran loops over j=2..22, but the exc/wid arrays are
    # only populated through state 21. With bounds checking enabled, the Fortran
    # crashes at state 22. Here I keep the driver loop intact and treat any
    # out-of-range state as contributing zero.

    exc = [
        0.0,
        0.0, 0.00444, 0.00765, 0.00964, 0.01084, 0.0137, 0.0151,
        0.0161, 0.0183, 0.020, 0.0230, 0.0315, 0.042,
        0.0151, 0.0161, 0.0166, 0.0181, 0.0193, 0.0206, 0.0235, 0.0315,
    ]
    wid = [
        0.0,
        0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00125,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00475, 0.009, 0.012,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00035, 0.00015, 0.004, 0.009,
    ]

    if state < 1 or state > 21:
        return 0.0, 0.0

    q2 = max(q2p, 0.0)
    qv2 = q2 + nu * nu
    # smwid = 0.0035
    smwid = 0.050
    # smwid = 0.001
    width = math.sqrt(smwid * smwid + wid[state] * wid[state])
    norm = width * math.sqrt(pi)

    nuel = q2 / (2.0 * (0.931494 * a))
    x = q2 / (2.0 * mp * nu)
    _ = x

    q2f = q2 / (0.1975 * 0.1975)
    if state == 18:
        exc[state] = min(0.0194 + 0.00016 * math.sqrt(q2f), 0.01955)

    nuex = nuel + exc[state]
    fs = _exp(-1.0 * (nu - nuex) ** 2 / width / width)
    fs = fs / norm
    if ((nu - nuex) / width) > 5.0:
        fs = 0.0

    ff = 0.0
    fft = 0.0
    if state <= 13:
        ff = nucffs12c(a, z, qv2, state)
    else:
        fft = nucffs12ct(a, z, qv2, state)

    w1 = 0.5 * (z * fft) ** 2 * fs
    wl = (z * ff) ** 2 * fs
    f1 = mp * w1
    fl = q2 * q2 / qv2 / nu * wl if qv2 > 0.0 and nu > 0.0 else 0.0
    return f1, fl


# -----------------------------------------------------------------------------
# csfitcomp.f
# -----------------------------------------------------------------------------


def csfitcomp(w2: float, q2: float, a: float, z: float, xvalc: Sequence[float], kind: int) -> Tuple[float, float]:
    """CSFITCOMP(w2,q2,A,Z,XVALC,type,sigt,sigL).

    Returns (sigt, sigl), corresponding to the Fortran outputs.
    """
    psimin = -2.3
    psimax = 5.0
    nbins = 220
    mp = 0.938272
    mp2 = mp * mp

    x = q2 / abs(w2 - mp2 + q2)
    dpsi = (psimax - psimin) / float(nbins)
    kappa2 = 1.0 + 4.0 * x * x * mp2 / q2

    int1 = 0.0
    int2 = 0.0
    for i in range(1, nbins + 1):
        psip = psimin + dpsi * (i - 1)
        fy1 = 1.5576 / (1.0 + 1.7720**2 * (psip + 0.3014) ** 2) / (1.0 + _exp(-2.4291 * psip))
        fy2 = xvalc[7] / (1.0 + xvalc[8] ** 2 * (psip + xvalc[9]) ** 2) / (1.0 + _exp(xvalc[10] * psip)) * (1.0 - abs(psip) / psimax) ** 2.0 * (1.0 + abs(psip) / psimax) ** 2.0
        if psip > psimax:
            fy2 = 0.0
        fy2 = max(0.0, fy2)
        int1 += fy1
        int2 += fy2
    _ = int1
    rat = 1.0 / int2 / dpsi

    f1i, f2i, fli = gsmearing(z, a, w2, q2, xvalc)
    if fli < 0.0:
        fli = 0.0

    r = fli / (2.0 * x * f1i) if x != 0.0 and f1i != 0.0 else 0.0
    _ = r

    f1qe, f2qe = qenuc21off(z, a, q2, w2, xvalc)
    f1qe *= rat
    f2qe *= rat
    flqe = kappa2 * f2qe - 2.0 * x * f1qe
    if flqe < 0.0:
        flqe = 0.0

    f1mec = mec2021(z, a, w2, q2, xvalc)
    flmec = 0.0
    f2mec = 2.0 * x * f1mec / kappa2

    if kind == 1:
        f1 = f1i + f1qe + f1mec
        f2 = f2i + f2qe + f2mec
        fl = fli + flqe
    elif kind == 2:
        f1 = f1qe
        f2 = f2qe
        fl = flqe
    elif kind == 3:
        f1 = f1i
        f2 = f2i
        fl = fli
    elif kind == 4:
        f1 = f1mec
        f2 = f2mec
        fl = flmec
    elif kind == 5:
        f1 = f1qe + f1mec
        f2 = f2qe + f2mec
        fl = flqe + flmec
    else:
        raise ValueError(f"Unsupported CSFITCOMP kind={kind}")

    if fl < 0.0:
        fl = 0.0

    sigt = f1
    sigl = fl / (2.0 * x)
    _ = f2
    return sigt, sigl




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
    """Calculate RL/RT response pieces for one kinematic point.

    Parameters
    ----------
    qv : float
        Three-momentum transfer |q| in GeV.
    nu : float
        Energy transfer in GeV.
    a : float
        Nuclear mass number A.
    z : float
        Nuclear charge Z.
    xvalc : sequence of float
        Fit parameter table used by the original Fortran code.

    Returns
    -------
    dict or None
        Dictionary with the response outputs. Returns None when q2 <= 0 or
        nu == 0, matching the original Fortran behavior.
    """
    q2 = qv * qv - nu * nu
    # if q2 <= 0.0 or nu == 0.0:
    if q2 < 0.0 or nu < 0.0:
        return None

    nuel = q2 / (2.0 * (0.931494 * a))
    ex = nu - nuel
    w2 = MP_MAIN * MP_MAIN + 2.0 * MP_MAIN * nu - q2
    xb = q2 / (2.0 * MP_MAIN * nu)

    # total, except nuclear states
    f1, fl = csfitcomp(w2, q2, a, z, xvalc, 1)
    fl = 2.0 * xb * fl
    rttot = 2.0 / MP_MAIN * f1 / 1000.0
    rltot = qv * qv / q2 / 2.0 / MP_MAIN / xb * fl / 1000.0

    # quasi elastic peak
    f1, fl = csfitcomp(w2, q2, a, z, xvalc, 2)
    fl = 2.0 * xb * fl
    rtqe = 2.0 / MP_MAIN * f1 / 1000.0
    rlqe = qv * qv / q2 / 2.0 / MP_MAIN / xb * fl / 1000.0

    # inelastic peak
    f1, fl = csfitcomp(w2, q2, a, z, xvalc, 3)
    fl = 2.0 * xb * fl
    rtie = 2.0 / MP_MAIN * f1 / 1000.0
    rlie = qv * qv / q2 / 2.0 / MP_MAIN / xb * fl / 1000.0

    # transverse enhancement
    f1, fl = csfitcomp(w2, q2, a, z, xvalc, 4)
    fl = 2.0 * xb * fl
    rte = 2.0 / MP_MAIN * f1 / 1000.0
    rle = 0.0
    _ = fl

    # nuclear states
    flns = 0.0
    f1ns = 0.0
    for state in range(2, 23):
        f1_state, fl_state = nuc12sf(z, a, nu, q2, state)
        flns += fl_state
        f1ns += f1_state
    rtns = 2.0 / MP_MAIN * f1ns / 1000.0
    rlns = qv * qv / q2 / 2.0 / MP_MAIN / xb * flns / 1000.0
    if rlns <= 1.0e-40:
        rlns = 0.0
    if rtns <= 1.0e-40:
        rtns = 0.0

    rltot += rlns
    rttot += rtns

    return {
        "qv": qv,
        "q2": q2,
        "ex": ex,
        "nu": nu,
        "rttot": rttot,
        "rltot": rltot,
        "rtqe": rtqe,
        "rlqe": rlqe,
        "rtie": rtie,
        "rlie": rlie,
        "rte": rte,
        "rle": rle,
        "rtns": rtns,
        "rlns": rlns
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
        # row = calculate_response_point(float(qv), float(nu), a=12, z=6, xvalc=xvalc)
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
    # rows = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    # for col in ["rttot", "rltot", "rtqe", "rlqe",
    #     "rtie", "rlie", "rte", "rle", "rtns", "rlns"]:
    #     rows[col] = rows[col] * z / 6 # scale from carbon to target nucleus
    # return rows


# -----------------------------------------------------------------------------
# Cross-section driver logic from qemodplot_norm.f / vcoul.f / nuccs12cs.f
# -----------------------------------------------------------------------------


def vcoul(a: float, z: float) -> float:
    """Return the effective Coulomb potential V used by the Fortran code.

    Parameters
    ----------
    a : float
        Nuclear mass number A.
    z : float
        Nuclear charge Z.

    Notes
    -----
    The original Fortran currently hard-codes the carbon value `V = 0.0031`
    GeV, even though the older comments describe a geometric estimate.
    This function follows that active Fortran behavior exactly.
    """
    _ = (a, z)
    # return 0.0031
    return 0.0074
    # if a == 40 and z == 20:
    #     return 0.0074
    # else:
    #     return 0.0031



def nuccs12cs(z: float, a: float, beam_energy: float, scattered_energy: float,
              scattering_angle_deg: float, state: int) -> float:
    """Return the 12C narrow-state cross section contribution in microbarn.

    Parameters
    ----------
    z : float
        Nuclear charge Z.
    a : float
        Nuclear mass number A.
    beam_energy : float
        Incoming electron energy E in GeV.
    scattered_energy : float
        Outgoing electron energy E' in GeV.
    scattering_angle_deg : float
        Electron scattering angle in degrees.
    state : int
        Narrow-state index. The Fortran driver sums over 2..21.
    """
    mp = 0.93827
    alpha = 7.29735e-03
    pi = 3.14159
    radcon = 0.0174533

    excitation_energy = [
        0.0, 0.00444, 0.00765, 0.00964, 0.01084, 0.0137, 0.0151,
        0.0161, 0.0183, 0.020, 0.0230, 0.0315, 0.042, 0.0151, 0.0161,
        0.0166, 0.0181, 0.0193, 0.0206, 0.0235, 0.0315, 0.01271,
    ]
    state_width = [
        0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00125,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00475, 0.009, 0.012,
        0.00002, 0.00002, 0.00002, 0.0002, 0.00035, 0.00015, 0.004,
        0.009, 0.0002,
    ]
    if not (1 <= state <= 22):
        return 0.0

    # smearing_width = 0.00048
    smearing_width = 0.002 # FIXME: smwidth
    total_width = math.sqrt(smearing_width * smearing_width + state_width[state - 1] * state_width[state - 1])
    normalization = total_width * math.sqrt(pi)

    sin_half_theta = math.sin(radcon * scattering_angle_deg / 2.0)
    sin2 = sin_half_theta * sin_half_theta
    cos2 = 1.0 - sin2
    if cos2 <= 0.0:
        return 0.0
    tan2 = sin2 / cos2

    energy_transfer = beam_energy - scattered_energy
    q2 = 4.0 * beam_energy * scattered_energy * sin2
    q2v = q2 + energy_transfer * energy_transfer
    if q2v <= 0.0 or q2 <= 0.0:
        return 0.0

    scattered_energy_elastic = a * 0.931494 * beam_energy / (a * 0.931494 + 2.0 * beam_energy * sin2)
    elastic_energy_transfer = beam_energy - scattered_energy_elastic

    state_excitation = excitation_energy[state - 1]
    q2f = q2 / 0.1975 / 0.1975
    if state == 18:
        state_excitation = min(0.0193 + 0.00015 * math.sqrt(max(q2f, 0.0)), 0.0195)

    excited_energy_transfer = elastic_energy_transfer + state_excitation
    q2_state = 4.0 * beam_energy * (beam_energy - excited_energy_transfer) * sin2
    q2v_state = q2_state + excited_energy_transfer * excited_energy_transfer
    if q2v_state <= 0.0:
        return 0.0

    gaussian_weight = math.exp(-1.0 * (energy_transfer - excited_energy_transfer) ** 2 / total_width / total_width)
    gaussian_weight = gaussian_weight / normalization
    if ((energy_transfer - excited_energy_transfer) / total_width) > 4.0:
        gaussian_weight = 0.0

    longitudinal_ff = 0.0
    transverse_ff = 0.0
    if state <= 13:
        longitudinal_ff = nucffs12c(a, z, q2v_state, state)
    else:
        transverse_ff = nucffs12ct(a, z, q2v_state, state)

    wl2 = (z * longitudinal_ff) ** 2
    wt2 = (z * transverse_ff) ** 2

    mott = 0.3894e3 * alpha ** 2 / 4.0
    if beam_energy <= 0.0 or sin2 <= 0.0:
        return 0.0
    mott = mott * cos2 / beam_energy / beam_energy / sin2 / sin2
    recoil = a * mp / 1.007276 / (12.0 * mp + beam_energy * (1.0 - math.cos(radcon * scattering_angle_deg)))

    narrow_state_cross_section = 1000.0 * mott * recoil * (
        q2 * q2 / q2v / q2v * wl2 + (q2 / 2.0 / q2v + tan2) * wt2
    )
    return gaussian_weight * narrow_state_cross_section



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

