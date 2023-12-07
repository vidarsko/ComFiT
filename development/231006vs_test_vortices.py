import comfit as cf
import numpy as np
import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201,gamma=0.08)
V_trap = bec.conf_harmonic_potential(100)
stirrer_radius = 20
stirrer_velocity = 0.6
t=0
#position = bec.update_circular_stirrer_position(stirrer_radius,t,stirrer_velocity)

#V_stirr = bec.gaussian_stirring_potential(4,3,position)

bec.V_ext = V_trap #+V_stirr

bec.conf_initial_condition_Thomas_Fermi()


#plt.figure(1)
#plt.imshow(np.abs(bec.psi))
#plt.colorbar()
bec.evolve_relax_BoseEinsteinCondensate(10)



#plt.imshow(np.abs(bec.psi))
# plt.show()

ax = plt.figure().add_subplot(111)

for n in range(20):
    print(n)
    bec.evolve_dGPE_with_stirrer(100,4,0.8,20,0.6)

    rho = bec.calc_vortex_density()

    plt.gcf().clear()
    # bec.plot_angle_field(np.angle(bec.psi))
    bec.plot_field(rho,colormap=cf.tool_colormap_bluewhitered(),cmin=-0.2,cmax=0.2)

    plt.draw()
    plt.pause(0.05)

