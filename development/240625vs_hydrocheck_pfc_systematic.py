import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import copy


# pfc = cf.PhaseFieldCrystal2DTriangular(30, 18)
pfc0 = cf.PhaseFieldCrystal2DSquare(30,30)
pfc0.dt=0.05

# This creates a standard orientation of the crystal
pfc0.conf_PFC_from_amplitudes()

# Create the rotated field
psi_rotated = pfc0.calc_PFC_from_amplitudes(rotation=np.pi/6)

# Specify the region centered at the mid position with radius 5 a0.
inclusion_region = pfc0.calc_region_disk(pfc0.rmid, 6*pfc0.a0)


# Set the rotated field in the inclusion region
pfc0.psi[inclusion_region] = psi_rotated[inclusion_region]
pfc0.psi_f = sp.fft.fftn(pfc0.psi)

fig = plt.figure(figsize=(15,5))
axs = fig.subplots(1,2)

tau = 10
pfc0.evolve_PFC(round(tau/pfc0.dt)) #Smooth the interface

delta_t = 10


for rho0exp, gammaSexp in zip(range(-6,3),range(-6,3)):
    pfc = copy.deepcopy(pfc0)
    method = 'ETD2RK'
    ID = f'{rho0exp}_{gammaSexp}_{method}'

    pfc.evolve_PFC_hydrodynamic(1, rho0=2**rho0exp, gamma_S=2**gammaSexp, method=method)

    for n in range(27):
        
        pfc.plot_field(pfc.psi[0], colorbar = False, ax=axs[0])

        f_f = pfc.calc_stress_divergence_f(pfc.psi_f[0])
        f = sp.fft.ifftn(f_f, axes =( range ( - pfc . dim , 0) ))
        f = np.real(f)/(pfc.el_mu/pfc.a0)

        pfc.plot_field(np.sqrt(f[0]**2 + f[1]**2), colorbar = True, ax=axs[1],
                        vlim=[0,7e-3])


        cf.tool_save_plot(n, ID=ID)
        fig.clear()     
        fig = plt.figure(figsize=(15,5))
        axs = fig.subplots(1,2)

        pfc.evolve_PFC_hydrodynamic(round(delta_t/pfc.dt), rho0=2**rho0exp, gamma_S=2**gammaSexp, method=method)

    cf.tool_make_animation_gif(n,name=f'hydrocheck_pfc_rho0_{rho0exp}_gammaS_{gammaSexp}_{method}.gif', ID=ID)