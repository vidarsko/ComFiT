import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

fig = plt.figure(figsize=(10, 10))

# Base System class instance
bs = cf.BaseSystem(2, xlim=[-3,3], ylim=[-3,3])
field = bs.x - bs.y
ax1 = fig.add_subplot(2, 2, 1) 
bs.plot_field(field, ax=ax1)

# # Quantum Mechanical System 
qm = cf.QuantumMechanics(3, xRes=41, yRes=41, zRes=41)
qm.conf_initial_condition_Gaussian(initial_velocity=[0,0.1,0.3])
qm.evolve_schrodinger(200)
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
qm.plot_complex_field(qm.psi, ax=ax2)

# Bose Einstein Condensate System
bec = cf.BoseEinsteinCondensate(3, xRes=41, yRes=41, zRes=41)
bec.conf_initial_condition_Thomas_Fermi()
bec.conf_insert_vortex_ring()
bec.evolve_relax(100)
vortex_nodes = bec.calc_vortex_nodes()
ax3 = fig.add_subplot(2, 2, 3, projection='3d')
bec.plot_field(abs(bec.psi), alpha = 0.2, ax=ax3, colorbar=False)
bec.plot_vortex_nodes(vortex_nodes,ax=ax3)

# Phase-field crystal system 
pfc = cf.PhaseFieldCrystal2DSquare(15,15)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)
dislocation_nodes = pfc.calc_dislocation_nodes()
ax4 = fig.add_subplot(2, 2, 4)
pfc.plot_field(pfc.psi,ax=ax4)
pfc.plot_dislocation_nodes(dislocation_nodes,ax=ax4,grid=False)

fig.tight_layout()

plt.show()
