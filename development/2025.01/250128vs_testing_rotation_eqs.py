
import comfit as cf
import numpy as np

nx=60

# Triangular lattice
ny = round(nx/np.sqrt(3))
pfc = cf.PhaseFieldCrystal2DTriangular(nx,ny)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
# eta = pfc.eta0
pfc.conf_PFC_from_amplitudes(eta, rotation=0.0*np.pi)
pfc.evolve_PFC(100)
# pfc.conf_create_polycrystal('four_grain')

dxxxpsi = pfc.ifft(pfc.dif[0] * pfc.dif[0] * pfc.dif[0] * pfc.fft(pfc.psi)).real
dxxypsi = pfc.ifft(pfc.dif[0] * pfc.dif[0] * pfc.dif[1] * pfc.fft(pfc.psi)).real
dxyypsi = pfc.ifft(pfc.dif[0] * pfc.dif[1] * pfc.dif[1] * pfc.fft(pfc.psi)).real
dyyypsi = pfc.ifft(pfc.dif[1] * pfc.dif[1] * pfc.dif[1] * pfc.fft(pfc.psi)).real

S111111 = pfc.ifft(pfc.fft(dxxxpsi**2)*pfc.calc_Gaussian_filter_f()).real
S111112 = pfc.ifft(pfc.fft(dxxxpsi*dxxypsi)*pfc.calc_Gaussian_filter_f()).real
Skkllmm = pfc.ifft(pfc.fft(dxxxpsi**2 + 3*dxxxpsi*dxxypsi + 3*dxxypsi**2 + dyyypsi**2 + dyyypsi**2)*pfc.calc_Gaussian_filter_f()).real

print("Mean of S111111:", np.mean(S111111))
print("Reference value of S111111:", pfc.S111111_ref)

op_x = (16/3*S111111 - 5/3*pfc.Skkllmm_ref)/(16/3*pfc.S111111_ref-5/3*pfc.Skkllmm_ref)
op_y = (16/3*S111112)/(16/3*pfc.S111111_ref-5/3*pfc.Skkllmm_ref)

op_unnormalized = op_x + 1j*op_y
op = np.tanh(np.abs(op_unnormalized))*np.exp(1j*np.angle(op_unnormalized))
op = op_unnormalized

op_f = pfc.fft(op)
dxtheta = np.imag(pfc.ifft(pfc.dif[0]*op_f)/op)
dytheta = np.imag(pfc.ifft(pfc.dif[1]*op_f)/op)


fig1 = pfc.plot_field(pfc.psi)
# fig2 = pfc.plot_angle_field(np.angle(op))
# fig2 = pfc.plot_field(np.angle(op)/6)
fig2 = pfc.plot_complex_field(op)


fig3 = pfc.plot_field(dxtheta)
fig4 = pfc.plot_field(dytheta)

fig = cf.tool_make_subplots(2,2,fig1,fig2,fig3,fig4)
fig.show()


