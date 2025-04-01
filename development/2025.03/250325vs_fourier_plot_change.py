import comfit as cf
import numpy as np

qm = cf.QuantumMechanics(1,xlim=[-10,10], xRes=100)

# def x_hat(psi):
#   return qm.x*psi

# def p_hat(psi):
#   return -1j*np.gradient(psi,qm.dx)

qm.conf_initial_condition_Gaussian(width = 0.1, position=0.0, initial_velocity=1)
# qm.psi = 0.5*np.sqrt(qm.calc_Gaussian(position=-1, width=1)) + 0.5*np.sqrt(qm.calc_Gaussian(position=3, width=1))

fig, axs = qm.plot_subplots(2,1)
qm.plot_complex_field(qm.psi,fig=fig,ax=axs[0])
qm.plot_complex_field(qm.psi_f,fig=fig,ax=axs[1],fourier=True)

qm.show(fig)