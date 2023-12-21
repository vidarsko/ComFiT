import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

pfc = cf.PhaseFieldCrystal2DTriangular(21,14, dt=0.1)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)


for n in range(20):
    pfc.evolve_PFC_hydrodynamic(100)

    # pfc.plot_angle_field(np.arctan2(pfc.psi[2],pfc.psi[1]),colorbar=False)
    pfc.plot_field(pfc.psi[2],colormap='winter')
    plt.draw()
    plt.pause(0.01)
