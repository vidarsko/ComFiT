import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


pfc = cf.PhaseFieldCrystal3DSimpleCubic(20,20,20)
pfc.plot_lib = 'plotly'
eta = pfc.calc_amplitudes_with_dislocation_ring()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(10)

Dnodes = pfc.calc_dislocation_nodes()
pfc.plot_dislocation_nodes(Dnodes)
plt.show()