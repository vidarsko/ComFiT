import comfit as cf
import numpy as np

bec = cf.BEC(2,xRes=51,yRes=51,gamma=0.05)
bec.conf_insert_vortex_dipole()
bec.evolve_relax_BEC(100)

for n in range(50):
    bec.evolve_dGPE(100)
    bec.plot_field(np.abs(bec.psi),colormap='winter',cmap_symmetric=False,
                   clims=[0,1])
    cf.tool_save_plot(n)

cf.tool_make_animation(n)