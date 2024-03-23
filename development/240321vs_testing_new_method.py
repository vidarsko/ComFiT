import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# pfc = cf.PhaseFieldCrystal2DSquare(20,20)
# pfc = cf.PhaseFieldCrystal1DPeriodic(20)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# pfc = cf.PhaseFieldCrystal2DTriangular(30,20,psi0=0.3)
# pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(1,1,1)
pfc = cf.PhaseFieldCrystal3DSimpleCubic(1,1,1)
print(pfc)
# pfc.conf_PFC_from_amplitudes()
# pfc.plot_field(pfc.psi)
# plt.show()

# bs = cf.BaseSystem(3,xlim=[0,4*np.pi*np.sqrt(2)],ylim=[0,4*np.pi*np.sqrt(2)],zlim=[0,4*np.pi*np.sqrt(2)])
# f = np.cos(bs.y/np.sqrt(2))*np.cos(bs.z/np.sqrt(2)) + np.cos(bs.x/np.sqrt(2))*np.cos(bs.z/np.sqrt(2)) + np.cos(bs.x/np.sqrt(2))*np.cos(bs.y/np.sqrt(2))
# bs.plot_field(f)
# plt.show()