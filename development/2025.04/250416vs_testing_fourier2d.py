import comfit as cf
import numpy as np
qm = cf.QuantumMechanics(2, xlim=[-10,15], ylim=[-15,10])
# qm = cf.QuantumMechanics(2)
qm.conf_initial_condition_Gaussian(position=[0,0], width=1, initial_velocity=[1,1])

fig, axs = qm.plot_subplots(1,2)
qm.plot_complex_field(qm.psi,  ax=axs[0])
qm.plot_complex_field(qm.psi_f,  ax=axs[1], fourier=True)
qm.show(fig)