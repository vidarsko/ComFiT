import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

bec = cf.BEC(2,xRes=101,yRes=101)

bec.conf_initial_condition_disordered()

# Create the figure and axes outside the loop
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

bec.plot_complex_field(bec.psi,ax)
for n in range(10):
    ax.cla()

    psi0 = bec.psi
    bec.evolve_dGPE(100)
    bec.plot_complex_field(bec.psi,ax)

    plt.pause(0.5)





