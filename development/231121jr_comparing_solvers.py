import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BEC(3,xRes=64,yRes=64,zRes=64,gamma=0.005,dt=0.01)

stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size =4
strength = .9

def V_t():
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * bec.t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * bec.t)
    pos_z = bec.zmid
    return  bec.calc_gaussian_stirring_potential(size, strength, [pos_x, pos_y,pos_z])

bec.V0 = V_t()


bec.set_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(20,'ETD2RK')



bec.set_time_dependent_potential(V_t)

bec.plot_field(np.abs(bec.psi)**2)
plt.show()

temp = np.copy(bec.psi)
temp_f = np.copy(bec.psi_f)

timesteps = int(200/bec.dt)

bec.evolve_dGPE(timesteps,'ETD2RK')

temp_1 = np.copy(bec.psi)
bec.psi = temp
bec.psi_f = temp_f

print("time is", bec.t)

bec.t = 0
bec.evolve_dGPE(timesteps,'ETD4RK')
print("time is", bec.t)

bec.plot_field(np.abs(temp_1 -bec.psi)**2)
plt.show()

print("Error:", np.sum(np.abs(temp_1-bec.psi)))