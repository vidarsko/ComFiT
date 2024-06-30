import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


# pfc = cf.PhaseFieldCrystal2DTriangular(30, 18)
pfc = cf.PhaseFieldCrystal2DSquare(30,30)
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(10,10,1)
pfc.dt=0.05

# This creates a standard orientation of the crystal
pfc.conf_create_polycrystal(type='four_grain')


fig = plt.figure(figsize=(15,5))
axs = fig.subplots(1,2)

delta_t = 10

rho0exp = -6
gammaSexp = -6
method = 'ETD2RK'
ID = f'{rho0exp}_{gammaSexp}_{method}'

pfc.evolve_PFC_hydrodynamic(1, rho0=2**rho0exp, gamma_S=2**gammaSexp, method=method)

strength = 0.01*pfc.el_mu/pfc.a0
pfc.external_force_density_f[0] = sp.fft.fftn(strength*np.sin(pfc.y/pfc.ymax*2*np.pi)+0*pfc.x)

for n in range(100):
    
    pfc.plot_field(pfc.psi[0], colorbar = False, ax=axs[0])

    f_f = pfc.calc_stress_divergence_f(pfc.psi_f[0])
    f = sp.fft.ifftn(f_f, axes =( range ( - pfc . dim , 0) ))
    f = np.real(f)/(pfc.el_mu/pfc.a0)

    pfc.plot_field(np.sqrt(f[0]**2 + f[1]**2), colorbar = True, ax=axs[1],
                    vlim=[0,0.02])


    cf.tool_save_plot(n, ID=ID)
    fig.clear()     
    fig = plt.figure(figsize=(15,5))
    axs = fig.subplots(1,2)

    pfc.evolve_PFC_hydrodynamic(round(delta_t/pfc.dt), rho0=2**rho0exp, gamma_S=2**gammaSexp, method=method)

cf.tool_make_animation_gif(n,name=f'hydrocheck_pfc_rho0_{rho0exp}_gammaS_{gammaSexp}_{method}.gif', ID=ID)