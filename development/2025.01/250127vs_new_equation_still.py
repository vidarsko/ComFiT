import comfit as cf
import numpy as np

pfc = cf.PhaseFieldCrystal2DSquare(40,40)

angle = 0.2*np.pi

eta = pfc.calc_amplitudes_with_dislocation_dipole()
# pfc.conf_PFC_from_amplitudes(eta, rotation=0.2)
pfc.conf_PFC_from_amplitudes(eta=eta,rotation=angle)
# pfc.conf_PFC_from_amplitudes(rotation=angle)
pfc.evolve_PFC(100)
psi = pfc.psi
print("Sii_eq from pfc:", pfc.Sii_eq)
Sii_eq = 4*pfc.A**2 + 8*pfc.B**2
print("Real Sii_eq:", Sii_eq)

print("Siijj_eq from pfc:", pfc.Siijj_eq)
S_iijj_eq = 4*pfc.A**2 + 16*pfc.B**2
print("Real Siijj_eq:", S_iijj_eq)


S1 = pfc.A**2 + 4*pfc.B**2
S2 = pfc.A**2 - 4*pfc.B**2

print("S1:", S1)
print("0.25 Siijj_eq", 0.25*S_iijj_eq)

print("S2:", S2)
print("Sii_eq - 0.75 Siijj_eq:", Sii_eq - 0.75*S_iijj_eq)

# q = np.array([[np.cos(angle), np.sin(angle)],
#               [-np.sin(angle), np.cos(angle)]])
# one_mode = np.exp( 1j*(q[0,0]*pfc.x + q[0,1]*pfc.y)) \
#     + np.exp(-1j*(q[0,0]*pfc.x + q[0,1]*pfc.y)) \
#     + np.exp( 1j*(q[1,0]*pfc.x + q[1,1]*pfc.y)) \
#     + np.exp(-1j*(q[1,0]*pfc.x + q[1,1]*pfc.y)) 
# Sii_eq_1m = 4

# p = np.zeros_like(q)
# p[0,:] = q[0,:] + q[1,:]
# p[1,:] = -q[0,:] + q[1,:]
# two_mode = np.exp( 1j*(p[0,0]*pfc.x + p[0,1]*pfc.y)) \
#          + np.exp(-1j*(p[0,0]*pfc.x + p[0,1]*pfc.y)) \
#          + np.exp( 1j*(p[1,0]*pfc.x + p[1,1]*pfc.y)) \
#          + np.exp(-1j*(p[1,0]*pfc.x + p[1,1]*pfc.y))
# Sii_eq_2m = 4*4

# psi = one_mode + two_mode
# Sii_eq = Sii_eq_1m + Sii_eq_2m
# S1 = 1 + 4
# S2 = 1 - 4

dxxpsi = pfc.ifft(pfc.dif[0]*pfc.dif[0]*pfc.fft(psi)).real
dxypsi = pfc.ifft(pfc.dif[0]*pfc.dif[1]*pfc.fft(psi)).real

# print("S1111_eq:", pfc.S1111_eq)
op_real = pfc.ifft(pfc.fft(2*dxxpsi*dxxpsi-3*S1)/S2*pfc.calc_Gaussian_filter_f()).real
op_imag = pfc.ifft(pfc.fft(2*dxxpsi*dxypsi)/S2*pfc.calc_Gaussian_filter_f()).real
op = op_real + 1j*op_imag

fig1 = pfc.plot_field(psi)
fig2 = pfc.plot_complex_field(op)

fig3 = pfc.plot_field(op_real)
fig4 = pfc.plot_field(op_imag)

fig = cf.tool_make_subplots(2,2,fig1,fig2,fig3,fig4)
fig.show()