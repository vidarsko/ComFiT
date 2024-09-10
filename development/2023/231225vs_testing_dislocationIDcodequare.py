import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint

dislocation_type=2
a = 0.3
b = 0.7
c = 0.1
pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(13, 13, 13)
eta = pfc.calc_amplitudes_with_dislocation_ring(
                    normal_vector=[a,b,c],
                    position=pfc.rmid,
                    radius=pfc.xmax/3,
                    dislocation_type=dislocation_type
                )
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(20)
#eta = pfc.calc_demodulate_PFC()
# pfc.plot_field(pfc.psi)
#pfc.plot_field(np.abs(eta[2]),cmap_symmetric=False)
# alpha = pfc.calc_dislocation_density()
Dnodes = pfc.calc_dislocation_nodes()
pprint(Dnodes)
pprint(pfc.a[dislocation_type-1])
# pfc.plot_field(alpha[0],colormap='bluewhitered')
# pfc.plot_field(pfc.psi,colorbar=False)

for dislocation in Dnodes:
    print((dislocation['Burgers_vector'] == pfc.a[dislocation_type-1]).all())

# pfc.plot_dislocation_nodes(Dnodes)
#print(eta[1])
# plt.show()