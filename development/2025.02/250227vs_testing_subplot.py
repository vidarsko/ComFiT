
import comfit as cf
import numpy as np

qm = cf.QuantumMechanics(2)
qm.conf_initial_condition_Gaussian(initial_velocity=[1,0])

fig, axs = qm.plot_subplots(2,3)

qm.plot_complex_field(qm.psi, fig=fig, ax=axs[0,0], title='Complex Field')
qm.plot_field(abs(qm.psi), fig=fig, ax=axs[0,1], title='Magnitude')
qm.plot_field(qm.psi.real, fig=fig, ax=axs[0,2], title='Real Part')
qm.plot_field(qm.psi.imag, fig=fig, ax=axs[1,0], title='Imaginary Part')
qm.plot_angle_field(np.angle(qm.psi), fig=fig, ax=axs[1,1], title='Phase')

# qm.show(fig)

qm.plot_save(1, fig)