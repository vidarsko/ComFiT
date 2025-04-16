import comfit as cf
import numpy as np
qm = cf.QuantumMechanics(1, xlim=[-10,10])
qm.conf_initial_condition_Gaussian(width=2, initial_velocity=5, position=0)

fig, axs = qm.plot_subplots(1,2)
qm.plot_complex_field(qm.psi,  ax=axs[0], title='Real Space')
qm.plot_complex_field(qm.psi_f, fourier=True,  ax=axs[1], title='Fourier Space')
qm.show(fig)