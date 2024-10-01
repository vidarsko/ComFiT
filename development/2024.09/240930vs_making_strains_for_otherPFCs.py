import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


# pfc = cf.PhaseFieldCrystal2DTriangular(30, 15, for_properties_calculation=False)
# pfc = cf.PhaseFieldCrystal2DSquare(2,2)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(1, 1, 1)
pfc = cf.PhaseFieldCrystal3DSimpleCubic(1, 1, 1)
# print(pfc.a)
# eta = pfc.calc_amplitudes_with_dislocation_dipole()
# pfc.conf_PFC_from_amplitudes(eta)
# pfc.evolve_PFC(100)
# pfc.plot_field(pfc.psi)
# stress =  pfc.calc_stress_tensor()
# fig = pfc.plot_field(stress[0])
# fig.show()
# print(pfc.calc_Gaussian_filter_f())
# print(pfc.a0)
# fig = pfc.plot_PFC()
# fig.show()
# print(np.sqrt(pfc.q[0,0]**2 + pfc.q[0,1]**2))
# print(pfc.q)