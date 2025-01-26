import comfit as cf
import numpy as np

nx=40
ny=round(nx/np.sqrt(3))
pfc = cf.PhaseFieldCrystal2DTriangular(nx,ny)
psi0 = pfc.calc_PFC_from_amplitudes()


N = 201
inner_product = np.zeros(N)
angles = np.linspace(0, np.pi/3, N)

for i in range(N):

    psi = pfc.calc_PFC_from_amplitudes(rotation = angles[i])
    
    dxpsi = pfc.ifft(pfc.dif[0]*pfc.fft(psi)).real
    dypsi = pfc.ifft(pfc.dif[1]*pfc.fft(psi)).real
    dxpsi0 = pfc.ifft(pfc.dif[0]*pfc.fft(psi0)).real
    dypsi0 = pfc.ifft(pfc.dif[1]*pfc.fft(psi0)).real

    inner_product[i] = pfc.calc_integrate_field(dxpsi*dxpsi0 + dypsi*dypsi0)/pfc.calc_integrate_field(psi*psi0)

    # inner_product[i] = 

import matplotlib.pyplot as plt

plt.figure()
plt.plot(angles, inner_product)
plt.xlabel('Angle (radians)')
plt.ylabel('Inner product')
plt.grid(True)
plt.show()


    

# pfc = cf.PhaseFieldCrystal2DSquare(40,40)
# eta = pfc.calc_amplitudes_with_dislocation_dipole()
# eta = pfc.calc_amplitudes_with_dislocation_dipole(x1=pfc.xmid, x2=pfc.xmid, y1 = pfc.ymax/3, y2 = 2*pfc.ymax/3, eta=eta, dislocation_type=1)
# pfc.conf_PFC_from_amplitudes(eta)
# # pfc.conf_PFC_from_amplitudes()
# pfc.conf_create_polycrystal(type='four_grain')
# pfc.evolve_PFC(100)

# fig = pfc.plot_field(pfc.psi)
# fig.show()

# S = pfc.calc_structure_tensor()
# S_f = pfc.fft(S)

# eta33 = 1/(6*pfc.A**2)*pfc.ifft(pfc.dif[0]**2*S_f[2] - 2*pfc.dif[0]*pfc.dif[1]*S_f[1] + pfc.dif[1]**2*S_f[0]).real
# eta33_1 = 0.003092104845449955/2
# print('Integral of eta33^2:', pfc.calc_integrate_field(eta33**2)/eta33_1)
# print('Integral of eta33:', pfc.calc_integrate_field(eta33))

# fig = pfc.plot_field(eta33)

# fig.show()

