import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

bec = cf.BEC(2,xRes=101,yRes=101)

bec.set_initial_condition_disordered()


for n in range(10):
    psi0 = bec.psi
    bec.evolve_dGPE(20)
    angle = np.angle(bec.psi)
    plt.clf()
    bec.plot_angle_field(angle)
    plt.pause(0.1)




