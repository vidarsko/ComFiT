import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
pfc.conf_PFC_from_amplitudes()

pfc.plot_field_in_plane(pfc.psi)
plt.show()