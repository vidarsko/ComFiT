import comfit as cf
import numpy as np
qm = cf.QuantumMechanics(3, xlim=[-10,15], ylim=[-15,10], zlim=[-11,12])
# qm = cf.QuantumMechanics(2)
qm.conf_initial_condition_Gaussian(position=[0,0,0], width=1, initial_velocity=[1,1,0])

fig, axs = qm.plot_subplots(1,2)
qm.plot_complex_field(qm.psi, fig=fig, ax=axs[0])
qm.plot_complex_field(qm.psi_f, fig=fig, ax=axs[1], fourier=True)
qm.show(fig)