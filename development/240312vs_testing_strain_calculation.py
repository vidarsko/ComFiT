import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

pfc = cf.PhaseFieldCrystal2DTriangular(45, 28)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)

pfc.evolve_PFC(100)
pfc.plot_field(pfc.psi)

strain = pfc.calc_strain_tensor()
strain = strain - np.mean(strain,axis=0)
pfc.plot_field(strain[2])

# structure_tensor = pfc.calc_structure_tensor()
# pfc.plot_field(structure_tensor[0])
plt.show()
