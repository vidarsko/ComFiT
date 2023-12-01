import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

qm = cf.QM(1,xRes=101,dx=1)

qm.set_initial_condition_gaussian(position=40,width=5)
#qm.set_harmonic_potential(0)

ymax = np.max(abs(qm.psi))**2
#qm.psi = qm.psi*np.exp(1j*qm.x/3)

for n in range(100):
    qm.evolve_schrodinger(100)
    qm.plot(ylim=ymax)
    # cf.tool_save_plot(n)
    plt.draw()
    plt.pause(0.01)
# cf.tool_make_animation(n)