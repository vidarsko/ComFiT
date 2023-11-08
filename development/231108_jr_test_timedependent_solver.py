import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BEC(2,xRes=64,yRes=64,gamma=0.01,dt=0.1)

stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size =4
strength = .9

def V_t(t):
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * t)
    return  bec.gaussian_stirring_potential(size, strength, [pos_x, pos_y])

bec.V_ext = V_t(0)
bec.set_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(20)

plt.imshow(np.abs(bec.psi)**2)
plt.show()

bec.evolve_time_dependent_ETDRK4(90000,V_t)

plt.imshow(np.abs(bec.psi)**2)
plt.colorbar()
plt.show()