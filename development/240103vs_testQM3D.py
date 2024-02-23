import os

import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf

import matplotlib.pyplot as plt


qm = cf.QuantumMechanics(3, xRes=71, yRes=71,zRes=71)

qm.conf_initial_condition_Gaussian(position=[0.5*qm.xmax,0.5*qm.ymax,0.35*qm.zmax],width=5)

qm.V_ext = qm.conf_harmonic_potential(0.001)

# qm.evolve_schrodinger(10)

# qm.plot_complex_field(qm.psi)

# plt.show()

# print(qm.psi)


for i in range(300):

    # qm.plot_field(abs(qm.psi))
    qm.plot_complex_field(qm.psi)

    cf.tool_save_plot(i)


    qm.evolve_schrodinger(10)

    # plt.draw()

    # plt.pause(0.01)
    


cf.tool_make_animation_movie(i)


