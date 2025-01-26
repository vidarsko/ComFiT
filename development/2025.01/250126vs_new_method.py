import comfit as cf
import numpy as np

pfc = cf.PhaseFieldCrystal2DSquare(40,40)

angle = 0.1*np.pi

eta = pfc.calc_amplitudes_with_dislocation_dipole()
# pfc.conf_PFC_from_amplitudes(eta, rotation=0.2)
# pfc.conf_PFC_from_amplitudes(eta=eta,rotation=angle)
pfc.conf_PFC_from_amplitudes(rotation=angle)
pfc.evolve_PFC(100)

width=0.5
psi = pfc.ifft(pfc.psi_f * np.exp(-(1-np.sqrt(pfc.calc_k2()))**2/(2*width**2))).real
psi = 4*psi/np.max(abs(psi))


dxpsi = pfc.ifft(pfc.dif[0]*pfc.fft(psi)).real
dypsi = pfc.ifft(pfc.dif[1]*pfc.fft(psi)).real

# pfc.plot_field(psi).show()

# print("S1111_eq:", pfc.S1111_eq)
op_real = pfc.ifft(pfc.fft(5 - 2/3*(dxpsi**4))*pfc.calc_Gaussian_filter_f()).real
# op_real = pfc.ifft(pfc.fft(1 + 2/3*(pfc.S1111_eq-dxpsi**4))*pfc.calc_Gaussian_filter_f()).real
op_imag = pfc.ifft(pfc.fft(-2/3*dxpsi**3*dypsi)*pfc.calc_Gaussian_filter_f()).real
op = op_real + 1j*op_imag


fig1 = pfc.plot_field(psi)
fig2 = pfc.plot_complex_field(op)

fig3 = pfc.plot_field(op_real)
fig4 = pfc.plot_field(op_imag)

fig = cf.tool_make_subplots(2,2,fig1,fig2,fig3,fig4)
fig.show()

# dxpsi = pfc.ifft(pfc.dif[0]*pfc.fft(psi)).real
# dypsi = pfc.ifft(pfc.dif[1]*pfc.fft(psi)).real
