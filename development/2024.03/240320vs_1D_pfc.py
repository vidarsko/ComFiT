import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

pfc = cf.PhaseFieldCrystal1DPeriodic(9,micro_resolution=[31])
pfc.conf_PFC_from_amplitudes()
pfc.conf_apply_strain(0.05)
for i in range(100):
    pfc.evolve_PFC(200)
    pfc.plot_field(pfc.psi)
    plt.pause(0.01)

plt.show()

