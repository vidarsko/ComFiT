import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


# pfc = cf.PhaseFieldCrystal2DTriangular(30, 18)
pfc = cf.PhaseFieldCrystal2DSquare(30,30)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(10,10,1)
pfc.dt=0.05

# This creates a standard orientation of the crystal
pfc.conf_create_polycrystal(type='four_grain')

orientation_field = pfc.calc_orientation_field()
pfc.plot_complex_field(orientation_field)

plt.show()
