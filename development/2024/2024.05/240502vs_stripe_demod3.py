import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

N=16
T=0.1

pfc = cf.PhaseFieldCrystal2DTriangular(60,23)
eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole(dipole_position=[0.33*pfc.xmax,pfc.ymax/2], dipole_vector=[pfc.xmax/4,0]))
theta=0.1
q0 = [np.cos(theta),np.sin(theta)]
psi1 = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))
theta=-0.1
q0 = [np.cos(theta),np.sin(theta)]
eta = eta*np.exp(1j*pfc.calc_angle_field_vortex_dipole(dipole_position=[0.66*pfc.xmax,pfc.ymax/2], dipole_vector=[-pfc.xmax/4,0]))
psi2 = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))
pfc.psi = psi1
pfc.psi[pfc.xRes//2:,:] = psi2[pfc.xRes//2:,:]
pfc.psi_f = np.fft.fft2(pfc.psi)
pfc.evolve_PFC(200)
# # pfc.plot_field(pfc.psi)
# # plt.show()

# pfc = cf.PhaseFieldCrystal2DTriangular(20,13)
# eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole())
# # q0 = [1,0]
# # q0=[np.sqrt(2)/2,np.sqrt(2)/2]
# # q0 = [0,1]
# q0 = [np.sqrt(3)/2,1/2]
# theta=np.pi/N/2
# q0 = [np.cos(theta),np.sin(theta)]
# pfc.psi = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))

eta0 = (np.max(pfc.psi)-np.min(pfc.psi))/2
print(eta0)



S = pfc.calc_structure_tensor()

fig = plt.figure()
axs = fig.subplots(3,3)

pfc.plot_field(S[0],ax=axs[0,0])
pfc.plot_field(S[1],ax=axs[0,1])
pfc.plot_field(S[2],ax=axs[0,2])

# def filter_field(psi):
#     psi_f = sp.fft.fftn(psi)
#     psi_f = psi_f*np.exp(-(pfc.calc_k2()-1))
#     return sp.fft.ifftn(psi_f)

# pfc.psi = filter_field(pfc.psi)
# pfc.psi_f = sp.fft.fftn(pfc.psi)

# fig = plt.figure()
# axs = fig.subplots(1,2)
# pfc.plot_field(pfc.psi,ax=axs[0])
# pfc.plot_field(filter_field(pfc.psi),ax=axs[1])
# plt.show()

def calc_B(pfc,strain):
    B = np.zeros((2,pfc.xRes,pfc.yRes))
    B[0] = sp.fft.ifftn(-pfc.dif[0]*sp.fft.fftn(strain[1]) + pfc.dif[1]*sp.fft.fftn(strain[0]))
    B[1] = sp.fft.ifftn(-pfc.dif[0]*sp.fft.fftn(strain[2]) + pfc.dif[1]*sp.fft.fftn(strain[1]))
    return np.real(B)

B = calc_B(pfc,S)/(2*eta0**2)

pfc.plot_field(B[0],ax=axs[1,0], vlim=[-0.005,0.005])
pfc.plot_field(B[1],ax=axs[1,1], vlim=[-0.005,0.005])

pfc.plot_field(pfc.psi,ax=axs[1,2])

Btot = np.sqrt(B[0]**2+B[1]**2)

pfc.plot_field(Btot,ax=axs[2,0])

plt.show()
