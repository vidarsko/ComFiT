import numpy as np
from scipy.fftpack import fftn, ifftn
from scipy.special import factorial
from collections import Counter

def tool_multinom(*args):
    M = len(args)
    if M < 2:
        return 1
    else:
        indices = np.array(args)
        counts = Counter(indices)
        prod = np.prod([factorial(count) for count in counts.values()])
        return 1 / prod

def levi_civita_symbol(i, j, k):
    if {i, j, k} != {0, 1, 2}:
        return 0  # Levi-Civita symbol is zero if indices are not 1, 2, 3

    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1  # Even permutation
    elif (i, j, k) in [(2, 1, 0), (0, 2, 1), (1, 0, 2)]:
        return -1  # Odd permutation
    else:
        return 0  # Repeated index

if __name__ == "__main__":
    print(tool_multinom(1,2,2))