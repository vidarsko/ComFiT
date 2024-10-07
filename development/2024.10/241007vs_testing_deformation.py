import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
import copy

for_properties_calculation = False

# pfc0 = cf.PhaseFieldCrystal2DTriangular(1,1, for_properties_calculation = for_properties_calculation)
# pfc0 = cf.PhaseFieldCrystal2DSquare(1,1, for_properties_calculation = for_properties_calculation)
pfc0 = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1, for_properties_calculation = for_properties_calculation)
# pfc0 = cf.PhaseFieldCrystal3DFaceCenteredCubic(1,1,1, for_properties_calculation = for_properties_calculation)
# pfc0 = cf.PhaseFieldCrystal3DSimpleCubic(1,1,1, for_properties_calculation = for_properties_calculation)

pfc0.conf_PFC_from_amplitudes()
pfc0.evolve_PFC(200)
f0 = pfc0.calc_free_energy()/pfc0.volume

strain_limit = 0.01
res = 101
strains = np.linspace(-strain_limit,strain_limit,res)

if pfc0.dim == 2:
    pure_shears = [[[0,a],[a,0]] for a in strains]
elif pfc0.dim == 3:
    pure_shears = [[[0,a,0],[0,a,0],[0,0,0]] for a in strains]

pure_shear_energies = np.zeros(res)
pure_shear_analytical = 2*pfc0.el_mu*(strains**2)

for i in range(res):
    pfc = copy.deepcopy(pfc0)
    pfc.conf_apply_distortion(pure_shears[i])
    f = pfc.calc_free_energy()/pfc.volume
    pure_shear_energies[i] = f-f0

if pfc0.dim == 2:
    pure_compressions = [[[a,0],[0,a]] for a in strains]
elif pfc0.dim == 3:
    pure_compressions = [[[a,0,0],[0,a,0],[0,0,0]] for a in strains]

pure_compression_energies = np.zeros(res)
pure_compression_analytical = (2*pfc0.el_lambda+2*pfc.el_mu+pfc0.el_gamma)*(strains**2)

for i in range(res):
    pfc = copy.deepcopy(pfc0)
    pfc.conf_apply_distortion(pure_compressions[i])
    f = pfc.calc_free_energy()/pfc.volume
    pure_compression_energies[i] = f-f0

if pfc0.dim == 2:
    volume_conserving_compressions = [[[a,0],[0,-a]] for a in strains]
elif pfc0.dim == 3:
    volume_conserving_compressions = [[[a,0,0],[0,-a,0],[0,0,0]] for a in strains]

volume_conserving_compression_energies = np.zeros(res)
volume_conserving_compression_analytical = (2*pfc0.el_mu + pfc0.el_gamma)*(strains**2)

for i in range(res):
    pfc = copy.deepcopy(pfc0)
    pfc.conf_apply_distortion(volume_conserving_compressions[i])
    f = pfc.calc_free_energy()/pfc.volume
    volume_conserving_compression_energies[i] = f-f0

fig, axs = plt.subplots(2, 3)

# Plot pure shear energies
# Plot pure shear energies
axs[0, 0].plot(strains, pure_shear_energies, label='Energies')
axs[0, 0].plot(strains, pure_shear_analytical, 'r--', label=r'$2 \mu \epsilon^2$')
axs[0, 0].set_title('Pure Shear')
axs[0, 0].set_xlabel('Strain')
axs[0, 0].set_ylabel('Energy Difference')
axs[0, 0].legend()
axs[0, 0].grid(True)

# Plot pure compression energies
axs[0, 1].plot(strains, pure_compression_energies, label='Energies')
axs[0, 1].plot(strains, pure_compression_analytical, 'r--', label=r'$(2 \lambda + 2 \mu + \gamma) \epsilon^2$')
axs[0, 1].set_title('Pure Compression')
axs[0, 1].set_xlabel('Strain')
axs[0, 1].set_ylabel('Energy Difference')
axs[0, 1].legend()
axs[0, 1].grid(True)

# Plot volume conserving compression energies
axs[0, 2].plot(strains, volume_conserving_compression_energies, label='Energies')
axs[0, 2].plot(strains, volume_conserving_compression_analytical, 'r--', label=r'$(2 \mu + \gamma) \epsilon^2$')
axs[0, 2].set_title('Volume Conserving Compression')
axs[0, 2].set_xlabel('Strain')
axs[0, 2].set_ylabel('Energy Difference')
axs[0, 2].legend()
axs[0, 2].grid(True)

# Plot differences between analytical and numerical solutions
axs[1, 0].semilogy(strains, np.abs(pure_shear_energies - pure_shear_analytical))
axs[1, 0].set_title('Difference in Pure Shear Energies')
axs[1, 0].set_xlabel('Strain')
axs[1, 0].set_ylabel('Absolute Energy Difference')
axs[1, 0].grid(True)

axs[1, 1].semilogy(strains, np.abs(pure_compression_energies - pure_compression_analytical))
axs[1, 1].set_title('Difference in Pure Compression Energies')
axs[1, 1].set_xlabel('Strain')
axs[1, 1].set_ylabel('Absolute Energy Difference')
axs[1, 1].grid(True)

axs[1, 2].semilogy(strains, np.abs(volume_conserving_compression_energies - volume_conserving_compression_analytical))
axs[1, 2].set_title('Difference in Volume Conserving Compression Energies')
axs[1, 2].set_xlabel('Strain')
axs[1, 2].set_ylabel('Absolute Energy Difference')
axs[1, 2].grid(True)


plt.tight_layout()
plt.show()
