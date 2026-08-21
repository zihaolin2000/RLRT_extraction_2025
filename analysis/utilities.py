
import numpy as np
from numpy.typing import ArrayLike
import numpy as np
from scipy.special import expit
from .presets import *

# ____________________utility functions____________________

def linear_model(x : ArrayLike, a : ArrayLike, b : ArrayLike) -> ArrayLike:
    return a * x + b

def v_coulomb(a : int = 12, z : int = 6) -> float:
    HBARC  = 0.197327      # in GeV.fm
    C_ASTE = 0.775
    r0     = 1.1*a**(1./3.) + 0.86*a**(-1./3.)

    # Coulomb potential at the center of the nucleus
    v0  = (3./2.)*ALPHA_FINE*HBARC*(z-1.)/r0  # in GeV

    # Average potential
    v  = C_ASTE*v0      # from Eur. Phys. J. A26 (2005) 167

    # # use experimentally determined value:
    # V = 0.0031 # Carbon
    # # V = 0.0081 # Iron

    # # V = 0.0

    return v

def pauli_blocking(ex : float) -> float:
    # ex in MeV
    if ex <=20:
        return np.exp(-73.3/ex)
    elif 20 < ex <= 140: 
        return 8.3714e-2 - 9.8343e-3 * ex + 4.1222e-4 * ex**2 - 3.4762e-6 * ex**3 + 9.3537e-9 * ex**4
    else:
        return np.exp(-24.2/ex)

def quasi_deuteron(ex : float, n : int = 6, a : int = 12, z : int = 6) -> float: # ex in GeV; return GD cross-section: GeV^-2
    ex = ex * 1e3 # convert ex to MeV
    if ex < 2.224:
        return 0
    else:    
        sigma = 397.8 * (n * z / a) * ((ex - 2.224)**(3/2)) * (ex**-3) * pauli_blocking(ex) # in mb
        return sigma * 0.1 * 0.1975**-2 # in GeV^-2

def dipole_form(q2s : ArrayLike) -> ArrayLike: # Q2 in GeV^2; return dipole form unitless
    # return 1 / ((1 + np.array(q2s) / 0.5)**5) # new value 2025 July 18
    return 1 / ((1 + np.asarray(q2s) / 0.5)**5) # new value 2025 July 18

def special_sigmoid(xs : ArrayLike, center : float = 0.12, width : float = 0.005) -> ArrayLike:
    z = (np.asarray(xs) - center) / width
    return expit(-z)

    # return np.where(xs <= center + 2 * width, 1 / (np.exp((np.array(xs) - center) / width) + 1), 0)

def rt_quasi_deuteron(nus : ArrayLike, q2s : ArrayLike, exs : ArrayLike) -> ArrayLike: # RT in MeV^-1, ex in GeV
    QDs=[]
    for ex in exs:
        QDs.append(quasi_deuteron(ex))
    QDs = np.asarray(QDs)
    GEs = dipole_form(q2s)
    GEs = np.asarray(GEs)
    # RTQD = GEs**2 * QDs * np.array(nus) / (2 * (np.pi**2) * ALPHA_FINE)
    RTQD = GEs**2 * QDs * np.asarray(nus) / (2 * (np.pi**2) * ALPHA_FINE)
    RTQD = RTQD * 1.5 * special_sigmoid(nus) # added 2025 Sep 23
    return RTQD*1e-3 # in MeV-1

def ratio_interpolated(x1, y1, x2, y2, eps=0.0, shrink=0.98):
    """
    Compute r = y2 / y1' on a shrunken shared x-range. The evaluation grid is taken
    from x2 values inside the shrunken overlap.

    Parameters
    ----------
    x1, y1 : array-like
    x2, y2 : array-like
    eps : float
        Small stabilizer added to denominator.
    shrink : float in (0,1]
        Fraction of the true overlap length to keep (centered). For example,
        0.98 keeps 98% of the overlap, trimming ~1% off each end.

    Returns
    -------
    r : np.ndarray
        Ratio y2 / interp(y1 at x_grid).
    x_grid : np.ndarray
        Monotonic array of x2 points inside the shrunken overlap.

    Raises
    ------
    ValueError if there is no overlap or no x2 points inside the shrunken range.
    """
    # to arrays & drop non-finite pairs
    x1, y1 = map(np.asarray, (x1, y1))
    x2, y2 = map(np.asarray, (x2, y2))
    m1 = np.isfinite(x1) & np.isfinite(y1)
    m2 = np.isfinite(x2) & np.isfinite(y2)
    x1, y1 = x1[m1], y1[m1]
    x2, y2 = x2[m2], y2[m2]

    if x1.size == 0 or x2.size == 0:
        raise ValueError("Empty x1/x2 after filtering non-finite values.")
    if not (0 < shrink <= 1):
        raise ValueError("shrink must be in (0, 1].")

    # sort x1 for interpolation
    o1 = np.argsort(x1)
    x1s, y1s = x1[o1], y1[o1]

    # true overlap
    lo = max(np.min(x1s), np.min(x2))
    hi = min(np.max(x1s), np.max(x2))
    if not (lo < hi):
        raise ValueError("No overlapping x-range between x1 and x2.")

    # shrink overlap symmetrically about its midpoint
    if shrink < 1.0:
        mid = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo) * shrink
        lo_s, hi_s = mid - half, mid + half
    else:
        lo_s, hi_s = lo, hi

    if not (lo_s < hi_s):
        raise ValueError("Shrunken overlap collapsed; try a larger 'shrink'.")

    # choose x-grid from x2 within shrunken overlap (sorted)
    mask = (x2 >= lo_s) & (x2 <= hi_s)
    if not np.any(mask):
        raise ValueError("No x2 points inside the shrunken overlap.")
    order2 = np.argsort(x2[mask])
    x_grid = x2[mask][order2]
    y2_grid = y2[mask][order2]

    # interpolate y1 -> x_grid (linear)
    y1_prime = np.interp(x_grid, x1s, y1s)

    # safe ratio
    denom = y1_prime + eps
    r = np.divide(y2_grid, denom,
                  out=np.full_like(y2_grid, np.nan, dtype=float),
                  where=np.isfinite(denom) & (denom != 0))

    return r, x_grid

def evaluate_norm_error(xs : ArrayLike, xs_err : ArrayLike, model : ArrayLike):
    ratio = model / xs
    ratio_err = ratio * np.sqrt(xs_err**2 / xs**2)
    ratio_weight = 1/ratio_err**2
    norm = np.sum(ratio * ratio_weight)/np.sum(ratio_weight)
    norm_err = np.sqrt(1/np.sum(ratio_weight))
    return norm, norm_err

