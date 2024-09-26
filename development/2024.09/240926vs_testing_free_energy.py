import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from scipy.integrate import simps

pfc = cf.PhaseFieldCrystal2DTriangular(1,1)
# pfc = cf.PhaseFieldCrystal2DSquare(1,1, micro_resolution=[7,7])
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
# pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(1,1,1)
# pfc = cf.PhaseFieldCrystal3DSimpleCubic(1,1,1)
pfc.conf_PFC_from_amplitudes()

N=300
t = np.linspace(0,N,N-1)
free_energy = np.zeros(len(t))
X,Y = np.meshgrid(pfc.x,pfc.y, indexing='ij')
for n in range(len(t)):
    pfc.evolve_PFC(1)
    free_energy_density, _ = pfc.calc_PFC_free_energy_density_and_chemical_potential()
    # print(free_energy_density.shape)
    # print(pfc.x.shape)
    # print(pfc.y.shape)
    # free_energy[n] = simps(simps(free_energy_density, pfc.y.flatten(), axis=1), pfc.x.flatten())
    free_energy[n] = pfc.calc_integrate_field(free_energy_density)

plt.plot(t,free_energy)
plt.show()

