import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

N=16
T=0.1

pfc = cf.PhaseFieldCrystal2DTriangular(40,13)
eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole(dipole_position=[0.33*pfc.xmax,pfc.ymax/2], dipole_vector=[pfc.xmax/4,0]))
theta=np.pi/N/2
q0 = [np.cos(theta),np.sin(theta)]
psi1 = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))
theta=-np.pi/N/2
q0 = [np.cos(theta),np.sin(theta)]
eta = eta*np.exp(1j*pfc.calc_angle_field_vortex_dipole(dipole_position=[0.66*pfc.xmax,pfc.ymax/2], dipole_vector=[-pfc.xmax/4,0]))
psi2 = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))
pfc.psi = psi1
pfc.psi[pfc.xRes//2:,:] = psi2[pfc.xRes//2:,:]
pfc.psi_f = np.fft.fft2(pfc.psi)
pfc.evolve_PFC(500)
# pfc.plot_field(pfc.psi)
# plt.show()



eta0 = (np.max(pfc.psi)-np.min(pfc.psi))/2
print(eta0)

fig = plt.figure()
axs = fig.subplots(2,3)
pfc.plot_field(pfc.psi, ax=axs[0,0])


def demodulate(pfc,q):
    return sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())


def orientation(pfc,T):
    qs = [[np.cos(theta),np.sin(theta)] for theta in np.linspace(0,np.pi-(np.pi)/N,N)]
    Z = 0
    p = np.zeros((len(qs),pfc.xRes,pfc.yRes))
    eta = np.zeros((len(qs),pfc.xRes,pfc.yRes),dtype=complex)
    for q,n in zip(qs,range(len(qs))):
        eta[n] = sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())
        p[n] = np.exp(abs(eta[n])/T)
        Z += p[n]
    
    Q = np.zeros((3,pfc.xRes,pfc.yRes))

    for q,n in zip(qs,range(len(qs))):
        prob = p[n]/Z
        Q[0] += prob*q[0]*q[0]
        Q[1] += prob*q[0]*q[1]
        Q[2] += prob*q[1]*q[1]
    
    return Q

Q = orientation(pfc,0.1)


def calc_strain(pfc,Q):
    strain = np.zeros((3,pfc.xRes,pfc.yRes))
    S = pfc.calc_structure_tensor()
    strain[0] = 1/2*Q[0] - 1/(4*eta0**2)*S[0]
    strain[1] = 1/2*Q[1] - 1/(4*eta0**2)*S[1]
    strain[2] = 1/2*Q[2] - 1/(4*eta0**2)*S[2]
    return strain 

strain = calc_strain(pfc,Q)

def calc_B(pfc,strain):
    B = np.zeros((2,pfc.xRes,pfc.yRes))
    B[0] = sp.fft.ifftn(-pfc.dif[0]*sp.fft.fftn(strain[1]) + pfc.dif[1]*sp.fft.fftn(strain[0]))
    B[1] = sp.fft.ifftn(-pfc.dif[0]*sp.fft.fftn(strain[2]) + pfc.dif[1]*sp.fft.fftn(strain[1]))
    return np.real(B)


pfc.plot_field(strain[0], ax=axs[0,1], title='strain[x,x]')
pfc.plot_field(strain[1], ax=axs[0,2], title='strain[x,y]')
pfc.plot_field(strain[2], ax=axs[1,0], title='strain[y,y]')

B = calc_B(pfc,strain)

pfc.plot_field(B[0], ax=axs[1,1], title='Bx')
pfc.plot_field(B[1], ax=axs[1,2], title='By')
plt.show()

    
