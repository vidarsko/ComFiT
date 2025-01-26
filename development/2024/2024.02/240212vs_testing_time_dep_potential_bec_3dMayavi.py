import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

# bec = cf.BoseEinsteinCondensate(3,xRes=101,yRes=101,zRes=101,gamma=0.02)

# bec.conf_initial_condition_Thomas_Fermi()

# bec.conf_external_potential(0)

# def V(t):
#     radius = 25*bec.a0
#     omega = 0.1
#     position = (bec.xmid+radius*np.cos(omega*t), bec.ymid+radius*np.sin(omega*t), bec.zmid)
#     return bec.calc_Gaussian(position=position, top = 2, width=10*bec.a0)

# # bec.evolve_relax(100)
# bec.conf_external_potential(V)

# for n in range(500):
#     mlab.clf()
#     bec.evolve_dGPE(20)

#     bec.plot_complex_field(bec.psi, plotting_lib='mayavi')
# # bec.plot_complex_field(bec.psi, plotting_lib='matplotlib')
#     # plt.draw()
#     # plt.pause(0.01)
#     # mlab.show()
#     mlab.savefig(f'plot_{n}.png')
#     mlab.close()
#     # cf.tool_save_plot_matplotlib(n)
cf.tool_make_animation_gif(211)
# plt.show()



# mlab.show()
