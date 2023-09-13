import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=101,yRes=101)
V_trap = bec.set_harmonic_potential(50)
stirrer_radius = 20
stirrer_velocity = 0.6
t=0
position = bec.position_update_circular_stirrer_2D(stirrer_radius,t,stirrer_velocity)

V_stirr = bec.gaussian_stiring_potential(4,3,position)

bec.V_ext = V_trap +V_stirr

bec.set_initial_condition_Thomas_Fermi()


#plt.figure(1)
#plt.imshow(np.abs(bec.psi))
#plt.colorbar()
bec.evolve_relax_BEC(10)



plt.imshow(np.abs(bec.psi))
# plt.show()

for i in range(1000):
    t += bec.dt
    position = bec.position_update_circular_stirrer_2D(stirrer_radius, t, stirrer_velocity)

    V_stirr = bec.gaussian_stiring_potential(4, 3, position)

    bec.V_ext = V_trap + V_stirr
    bec.evolve_dGPE(1)
    if i % 100 == 0:
        plt.imshow(np.abs(bec.psi))
        plt.draw()
        plt.pause(0.01)


