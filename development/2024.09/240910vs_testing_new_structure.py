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
pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# pfc.conf_PFC_from_amplitudes()
complex_field = 0.5+1j*np.sin(0.2*pfc.x)
fig = pfc.plot_complex_field(complex_field)
fig.show()
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

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)

# vector_field = np.zeros((2, pfc.xRes))
# vector_field[0] = np.sin(pfc.x)
# vector_field[1] = np.cos(pfc.x)

# fig = pfc.plot_vector_field(vector_field)
# fig.show()

## Field in plane

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# field = np.sin(pfc.x)

# pfc = cf.PhaseFieldCrystal2DTriangular(10,5)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# field = np.sin(pfc.x/pfc.a0)*np.sin(pfc.y/pfc.a0)*np.sin(pfc.z/pfc.a0)


# fig, ax = cf.plot_field_in_plane_matplotlib(pfc, field)
# plt.show()

# pfc.plot_field_in_plane(field)

## Complex field in plane

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# complex_field = 0.5+1j*np.sin(0.2*pfc.x)

# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# complex_field = 0.5+1j*np.sin(0.2*pfc.x)*np.sin(0.2*pfc.y)*np.sin(0.2*pfc.z)

# fig = pfc.plot_complex_field_in_plane(complex_field)
# fig, ax = cf.plot_complex_field_in_plane_matplotlib(pfc, complex_field)
# plt.show()

## Angle field in plane

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# angle_field = 0.5*pfc.x

# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# angle_field = 0.5*pfc.x*np.sin(0.2*pfc.y)*np.sin(0.2*pfc.z)

# cf.plot_angle_field_in_plane_matplotlib(pfc, angle_field)
# plt.show()
# pfc.plot_angle_field_in_plane(angle_field)



## Vector field in plane

# pfc = cf.PhaseFieldCrystal1DPeriodic(10)
# vector_field = np.zeros((2, pfc.xRes))

# pfc.plot_vector_field_in_plane(vector_field)

# bs = cf.BaseSystem(3,xRes=11,yRes=11,zRes=11)
# vector_field = np.array([bs.z+bs.x*np.cos(bs.y/5), bs.z+bs.y*np.sin(bs.x/5), -bs.z+bs.x*np.cos(bs.y/5)])


# bs = cf.BaseSystem(3,xRes=11,yRes=11,zRes=11)
# vector_field = np.array([bs.z+bs.x*np.cos(bs.y/5), bs.z+bs.y*np.sin(bs.x/5)])

# fig, ax = cf.plot_vector_field_in_plane_matplotlib(bs, vector_field, normal_vector=[0,1,1],position=[2,3,3])
# plt.show()
