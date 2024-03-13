import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

pfc = cf.PhaseFieldCrystal2DSquare(40, 40)

# This creates a standard orientation of the crystal
pfc.conf_PFC_from_amplitudes()


psi = pfc.calc_PFC_from_amplitudes(rotation=0.1)
pfc.plot_field(psi)
plt.show()