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

if __name__ == "__main__":
    print(tool_multinom(1,2,2))