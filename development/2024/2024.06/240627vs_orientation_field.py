import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


# pfc = cf.PhaseFieldCrystal2DTriangular(30, 18)
pfc = cf.PhaseFieldCrystal2DSquare(30,30)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(10,10,1)
# pfc.dt=0.05

# This creates a standard orientation of the crystal
pfc.conf_create_polycrystal(type='circular')
# pfc.conf_create_polycrystal(type='four_grain')
# pfc.evolve_PFC_hydrodynamic(1)
# pfc.plot_field(pfc.psi)


angle = np.zeros((pfc.xRes, pfc.yRes))
# complex_field = np.exp(1j * angle)
# pfc.plot_complex_field(complex_field)
# pfc.plot_angle_field(angle)
orientation_field = pfc.calc_orientation_field()
# pfc.plot_vector_field(orientation_field)
pfc.plot_orientation_field(orientation_field)
# angle = np.arctan2(orientation_field[1], orientation_field[0])
# pfc.plot_field(angle)

plt.show()
