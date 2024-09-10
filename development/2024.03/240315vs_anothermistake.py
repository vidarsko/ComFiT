import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


pfc = cf.PhaseFieldCrystal2DTriangular(34,20)
print(pfc.eta0)
# print(2*np.ones([2]+pfc.dims))hopw

eta = pfc.calc_amplitudes_with_dislocation_dipole()
print(eta)
pfc.conf_PFC_from_amplitudes(eta)

# pfc.plot_field(np.angle(eta[2]))
pfc.plot_complex_field(eta[0])
# pfc.plot_field(pfc.psi)
plt.show()

