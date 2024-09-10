import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BoseEinsteinCondensate(3,xRes=64,yRes=64,zRes=64,gamma=0.005,dt=0.1)

stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size =4
strength = .9

def V_t(t):
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * t)
    pos_z = bec.zmid
    return  bec.calc_gaussian_stirring_potential(size, strength, [pos_x, pos_y,pos_z])

#bec.conf_time_dependent_potential(V_t)

bec.V_ext = V_t(0)
bec.conf_initial_condition_Thomas_Fermi()
bec.evolve_relax(20)

bec.plot_field(np.abs(bec.psi)**2)
plt.show()

bec.evolve_time_dependent_ETDRK4(2000,V_t)

bec.plot_field(np.abs(bec.psi)**2)
plt.show()