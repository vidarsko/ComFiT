import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
pfc.evolve_PFC(500)
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
    delta_q = np.zeros((2,pfc.xRes,pfc.yRes))

    for q,n in zip(qs,range(len(qs))):
        prob = p[n]/Z
        eta_f = sp.fft.fftn(eta[n])
        # delta_q[0] = np.imag(sp.fft.ifftn(pfc.calc_Gaussian_filter_f()*sp.fft.fftn(sp.fft.ifftn(pfc.dif[0]*eta_f)/eta[n])))
        # delta_q[1] = np.imag(sp.fft.ifftn(pfc.calc_Gaussian_filter_f()*sp.fft.fftn(sp.fft.ifftn(pfc.dif[1]*eta_f)/eta[n])))
        delta_q[0] = np.imag(sp.fft.ifftn(pfc.dif[0]*eta_f)/eta[n])
        delta_q[1] = np.imag(sp.fft.ifftn(pfc.dif[1]*eta_f)/eta[n])

        Q[0] += abs(eta[n])**2*(q[0]+delta_q[0])*(q[0]+delta_q[0])
        Q[1] += abs(eta[n])**2*(q[0]+delta_q[0])*(q[1]+delta_q[1])
        Q[2] += abs(eta[n])**2*(q[1]+delta_q[1])*(q[1]+delta_q[1])

        # fig = plt.figure()
        # axs = fig.subplots(2,2)
        # pfc.plot_field(delta_q[0], ax=axs[0,0], title='delta_q[0]')
        # pfc.plot_field(delta_q[1], ax=axs[0,1], title='delta_q[1]')
        # pfc.plot_field(abs(eta[n]), ax=axs[1,0], title='abs(eta)')
        # pfc.plot_field(prob, ax=axs[1,1], title='prob')
        # plt.show()



    
    return Q

Q = orientation(pfc,0.001)

fig = plt.figure()
axs = fig.subplots(2,3)

pfc.plot_field(pfc.psi, ax=axs[0,0], title='psi')
pfc.plot_field(Q[0], ax=axs[0,1], title='Q[0]')
pfc.plot_field(Q[1], ax=axs[1,0], title='Q[1]')
pfc.plot_field(Q[2], ax=axs[1,1], title='Q[2]')

pfc.plot_field(Q[0]+Q[2], ax=axs[1,2], title='Q[0]+Q[2]')

Qeig = np.zeros((pfc.xRes,pfc.yRes,2,2))
Qeig[:,:,0,0] = Q[0]
Qeig[:,:,1,1] = Q[2]
Qeig[:,:,0,1] = Q[1]
Qeig[:,:,1,0] = Q[1]

eigv, q0 = np.linalg.eigh(Qeig)


psi_reconstructed = eta0*np.real(np.exp(1j*(q0[:,:,1,0]*pfc.x + q0[:,:,1,1]*pfc.y)))

pfc.plot_field(psi_reconstructed, ax=axs[0,2], title='psi_reconstructed')
plt.show()
