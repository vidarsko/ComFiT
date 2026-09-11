from typing import Tuple

# General packages
import numpy as np


def tool_generate_numbers_between(
        cmin: float,
        cmax: float
        ) -> Tuple[np.ndarray, float]:
    """Generate evenly spaced numbers between given minimum and maximum values.

    This function generates a sequence of numbers between cmin and cmax with
    appropriate step sizes based on the range magnitude. The step size is
    adjusted to maintain a readable number of ticks.

    Parameters
    ----------
    cmin : float
        Minimum value of the range.
    cmax : float
        Maximum value of the range.

    Returns
    -------
    tuple[np.ndarray, float]
        A tuple containing:
        - np.ndarray: Array of evenly spaced numbers between cmin and cmax
        - float: The exponent of the range magnitude (delta_exp)

    Raises
    ------
    ValueError
        If cmin is greater than cmax, or if either is NaN.

    Notes
    -----
    The step size is determined based on the range magnitude:
    - If range <= 5 units: step = 10^delta_exp
    - If range <= 10 units: step = 2×10^delta_exp
    - If range <= 15 units: step = 3×10^delta_exp
    - If range <= 20 units: step = 4×10^delta_exp
    """
    if np.isnan(cmin) or np.isnan(cmax):
        raise ValueError("cmin and cmax must not be NaN")

    if cmin > cmax:
        raise ValueError("cmin must be less than cmax")

    delta = (cmax - cmin)

    if delta == 0:
        return [0], 0, False

    delta_exp = np.floor(np.log10(max(abs(cmin), abs(cmax))))

    tick_min = np.ceil(cmin/10**delta_exp)*10**delta_exp
    tick_max = np.floor(cmax/10**delta_exp)*10**delta_exp

    number_of_steps = np.floor(delta/10**delta_exp)

    reiterated = False
    if number_of_steps < 1:
        cmid = (cmin + cmax) / 2
        numbers_between, delta_exp, reiterated = tool_generate_numbers_between(cmin-cmid, cmax-cmid)
        reiterated = True

    elif number_of_steps <= 5:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 10**delta_exp)
    elif number_of_steps <= 10:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 2*10**delta_exp)
    elif number_of_steps <= 15:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 3*10**delta_exp)
    elif number_of_steps <= 20:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 4*10**delta_exp)

    return numbers_between, delta_exp, reiterated
