import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(13, 13, 13)
eta = pfc.calc_amplitudes_with_dislocation_ring()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(10)

stress = pfc.calc_stress_tensor()

fig = pfc.plot_field_in_plane(stress[0], normal_vector=[1,1,1])
pfc.plot_field_in_plane(stress[0], normal_vector=[1,1,-1], fig=fig, colorbar=False)
fig.show()