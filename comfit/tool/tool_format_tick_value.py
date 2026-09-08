# General packages
import numpy as np


def tool_format_tick_value(
    val: float
    ) -> str:
    """Format a numeric value with appropriate SI prefix and significant figures.
    This function formats numeric values using SI prefixes and adjusts decimal places
    based on the magnitude of the number. For values that don't align with standard
    SI prefixes, scientific notation is used as a fallback.

    Parameters
    ----------
    val : float
        The numeric value to format.

    Returns
    -------
    str
        A string representation of the value with appropriate SI prefix and formatting.
    """
    if val == 0:
        return '0'

    # prefixes = {exp: f'×10{"".join(superscripts[digit] for digit in str(exp))}' for exp in range(-18, 19, 3)}
    # prefix_letters = {
    #      -15: 'f', -12: 'p', -9: 'n', -6: 'µ', -3: 'm',
    #     0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T'
    # }
    # for key in prefix_letters.keys():
    #     prefixes[key] = prefix_letters[key]


    abs_val = abs(val)
    exp = int(np.floor(np.log10(abs_val) / 3) * 3)  # Round to nearest power of 1000

    prefixes = {}

    if exp in prefixes:
        scaled_val = val / 10**exp
        # Format based on the magnitude of scaled_val
        if abs(scaled_val) < 10:
            return f'{scaled_val:.2f}{prefixes[exp]}'  # 2.00m
        elif abs(scaled_val) < 100:
            return f'{scaled_val:.1f}{prefixes[exp]}'  # 20.0m
        else:
            return f'{scaled_val:.0f}{prefixes[exp]}'  # 200m
    else:
        return f'{val:.0e}'  # Fallback to scientific notation
