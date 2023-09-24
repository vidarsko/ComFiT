import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

bec = cf.BEC(2,xRes=101,yRes=101)

bec.set_initial_condition_disordered()
bec.dt=0.1

# Create the figure and axes outside the loop
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111,projection='3d')

# bec.plot_complex_field(bec.psi,ax)
# bec.plot_field(bec.calc_k2(),ax)
# integration_factors_f = bec.calc_evolution_integrating_factors_dGPE_f()
# print(integration_factors_f[2][0,0])
# bec.plot_complex_field(integration_factors_f[2],ax)
# plt.show()
# plt.pause(10)


# plt.pause(5)
bec.plot_complex_field(bec.psi,ax)
for n in range(100):
    print(n)
    ax.cla()

    psi0 = bec.psi
    # bec.evolve_relax_BEC(100)
    bec.evolve_dGPE(100)
    bec.plot_complex_field(bec.psi,ax)
    # bec.plot_angle_field(np.angle(bec.psi))

    plt.pause(0.2)





