import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import mayavi.mlab as mlab
import numpy as np

print(mlab.get_engine().scenes)



bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31,gamma=0.01)

bec.conf_initial_condition_Thomas_Fermi()

bec.conf_external_potential(0)

def V(t):
    radius = 5*bec.a0
    omega = 0.1
    position = (bec.xmid+radius*np.cos(omega*t), bec.ymid+radius*np.sin(omega*t),bec.zmid)
    return bec.calc_Gaussian(position=position, top = 2, width=2*bec.a0)

# bec.evolve_relax(100)
bec.conf_external_potential(V)

# fig = mlab.figure(bgcolor=(1, 1, 1))

for n in range(1000):
    bec.evolve_dGPE(10)
    bec.plot_complex_field(bec.psi, plotting_lib='mayavi')
    # plt.pause(0.01)
    # plt.draw()
    # mlab.show()
    # mlab.draw()
    # mlab.close()
    # cf.tool_save_plot(n)
# cf.tool_make_animation_gif(n)
