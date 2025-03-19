
import comfit as cf
import numpy as np

qm = cf.QuantumMechanics(3,xRes=21,yRes=21,zRes=21)
# qm.plot_lib='matplotlib'

qm.psi = qm.y*np.exp(2*np.pi*1j*qm.x/qm.xmax)

fig, ax = qm.plot_complex_field(qm.psi, plot_method='phase_blob')
qm.show(fig)
