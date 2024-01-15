import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt


# 240110vs_testing_PFCmecheq.py

pfc = cf.PhaseFieldCrystal2DTriangular(21,12)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

u = pfc.calc_displacement_field_to_equilibrium()
# pfc.plot_vector_field(u)
# plt.show()