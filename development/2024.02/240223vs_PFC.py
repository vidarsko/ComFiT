import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import mayavi.mlab as mlab

# pfc = cf.PhaseFieldCrystal2DTriangular(30,20)
pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(10,10,10)
# eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes()
pfc.plot_field(pfc.psi, plotting_lib='mayavi')
mlab.show()
# plt.show()