# General packages
import numpy as np
from scipy.fftpack import fftn, ifftn
from scipy.special import factorial
from collections import Counter

def tool_multinom(
        *args: int
        ) -> float:
    """Calculate the multinomial coefficient.

    Parameters
    ----------
    \*args : int
        Sequence of integers representing the counts of each category.

    Returns
    -------
    float
        The multinomial coefficient value.
    """
    M = len(args)
    if M < 2:
        return 1
    else:
        indices = np.array(args)
        counts = Counter(indices)
        prod = np.prod([factorial(count) for count in counts.values()])
        return 1 / prod

def levi_civita_symbol(i, j, k):
    """Calculate the Levi-Civita symbol (antisymmetric tensor).

    The Levi-Civita symbol is defined as:
        1 if (i,j,k) is an even permutation of (0,1,2)
        -1 if (i,j,k) is an odd permutation of (0,1,2)
        0 if any indices are repeated or not in {0,1,2}
    
    Parameters
    ----------
    i : int
        First index
    j : int
        Second index  
    k : int
        Third index

    Returns
    -------
    int
        Value of Levi-Civita symbol: 1, -1, or 0

    Examples
    --------
    >>> levi_civita_symbol(0,1,2)
    1
    >>> levi_civita_symbol(1,0,2) 
    -1
    >>> levi_civita_symbol(0,0,1)
    0
    """
    if {i, j, k} != {0, 1, 2}:
        return 0  # Levi-Civita symbol is zero if indices are not 0, 1, 2

    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1  # Even permutation
    elif (i, j, k) in [(2, 1, 0), (0, 2, 1), (1, 0, 2)]:
        return -1  # Odd permutation
    else:
        return 0  # Repeated index

if __name__ == "__main__":
    print(tool_multinom(1,2,2))