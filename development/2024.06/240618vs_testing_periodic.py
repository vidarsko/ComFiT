import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


bs = cf.BaseSystem(1,xlim=[-1,1], xRes=10)
print(bs.x)
print(bs.k[0])

print(bs.k[0]*bs.x)