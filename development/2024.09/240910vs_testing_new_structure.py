import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# pfc = cf.PhaseFieldCrystal2DTriangular(10,5)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# pfc.conf_PFC_from_amplitudes()
# complex_field = 0.5+1j*np.sin(0.2*pfc.x)
# fig = pfc.plot_complex_field(complex_field)
# fig.show()
# cf.plot_field_matplotlib(pfc, pfc.psi)

# cf.plot_complex_field_matplotlib(pfc, complex_field)
# plt.show()
# fig.show()
# plt.show()

## Angle plot

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# pfc = cf.PhaseFieldCrystal2DTriangular(10,5)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)


# angle_field = 0.5*pfc.x
# # pfc.plot_angle_field(angle_field)

# cf.plot_angle_field_matplotlib(pfc, angle_field)
# plt.show()


## Vector fields

pfc = cf.PhaseFieldCrystal1DPeriodic(10)

vector_field = np.zeros((2, pfc.xRes))
vector_field[0] = np.sin(pfc.x)
vector_field[1] = np.cos(pfc.x)

fig = pfc.plot_vector_field(vector_field)
fig.show()