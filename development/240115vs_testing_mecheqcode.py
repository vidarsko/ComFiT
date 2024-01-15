import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DTriangular(21,13)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

# pfc.plot_field(pfc.psi)
# plt.show()
# OK

u = pfc.calc