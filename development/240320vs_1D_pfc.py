import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

pfc = cf.PhaseFieldCrystal1DPeriodic(9,micro_resolution=[31])
pfc.conf_PFC_from_amplitudes()
pfc.evolve_PFC(100)
pfc.a0=1
pfc.x = pfc.x-3*2*np.pi
pfc.plot_field(pfc.psi, xlabel='x', ylabel= r'$\psi$',xlim = [-6*np.pi,6*np.pi])
plt.show()
