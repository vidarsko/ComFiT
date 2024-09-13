import comfit as cf

import numpy as np

import matplotlib.pyplot as plt


qm = cf.QuantumMechanics(1,xRes=101,dx=1)


qm.conf_initial_condition_Gaussian(position=30,width=5)

qm.V_ext = qm.conf_harmonic_potential(0.0001)



ymax = np.max(abs(qm.psi))**2


for n in range(1000):
    qm.plot(ylim=ymax)

    cf.tool_save_plot_matplotlib(n)

    #plt.draw()

    #plt.pause(0.01)

    qm.evolve_schrodinger(10)

    ymax = np.max([ymax, np.max(abs(qm.psi)) ** 2])

cf.tool_make_animation_movie(n)