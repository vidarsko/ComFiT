import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# pfc = cf.PhaseFieldCrystal2DSquare(20,20)
# pfc = cf.PhaseFieldCrystal1DPeriodic(20)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
pfc = cf.PhaseFieldCrystal2DTriangular(1,1,t=1)
print(pfc)