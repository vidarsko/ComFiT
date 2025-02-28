
import comfit as cf
import numpy as np

qm = cf.QuantumMechanics(2)
qm.conf_initial_condition_Gaussian(initial_velocity=[1,0])
qm.psi = 0.23*qm.psi

fig, ax = qm.plot_field(qm.psi.imag, title='Imaginary Part')

qm.show(fig)

# qm.plot_save(1, fig)