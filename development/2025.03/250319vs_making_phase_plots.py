
import comfit as cf
import numpy as np

qm = cf.QuantumMechanics(3,xRes=21,yRes=21,zRes=21)
# qm.plot_lib='matplotlib'

qm.psi = np.exp(2*np.pi*1j*qm.x/qm.xmax + 0.5*np.pi*1j*qm.y/qm.ymax + 0.25*np.pi*1j*qm.z/qm.zmax)


fig, ax = qm.plot_subplots(1,1)
ax = ax[0]
# qm.plot_complex_field(qm.psi, plot_method='phase_blob', fig=fig, ax=ax)
qm.plot_complex_field(qm.psi, plot_method='phase_angle', fig=fig, ax=ax)
# qm.plot_complex_field_in_plane(qm.psi, normal_vector=[0,0,1], fig=fig, ax=ax)
qm.show(fig)
