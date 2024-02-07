import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
from mayavi import mlab

# bs = cf.BaseSystem(3,xRes=51,dx=0.1,yRes=51,dy=0.1,zRes=51,dz=0.1)
# z = bs.x**2 + bs.y**2 + bs.z**2

pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(3,3,3)
pfc.conf_PFC_from_amplitudes()

# bs.plot_field(z, number_of_layers=1, plotting_lib='mayavi')
# pfc.plot_field(pfc.psi, plotting_lib='matplotlib')
# plt.show()
pfc.plot_field(pfc.psi, plotting_lib='mayavi', number_of_layers=2)
mlab.show()