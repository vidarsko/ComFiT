import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


pfc = cf.PhaseFieldCrystal3DSimpleCubic(3,3,3)

pfc.conf_PFC_from_amplitudes(rotation=[0,0,0.2])

pfc.plot_field(pfc.psi)

plt.show()