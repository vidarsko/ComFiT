import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

pfc = cf.PhaseFieldCrystal2DTriangular(20,13)
# eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole())
# pfc.psi = np.real(eta*np.exp(1j*pfc.x))

pfc.psi = np.random.rand(pfc.xRes,pfc.yRes)-0.5

pfc.psi_f = np.fft.fft2(pfc.psi)


pfc.evolve_PFC(500)

def demodulate(pfc,q):
    return sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())


def orientation(pfc,T):
    qs = [[1,0],[np.sqrt(3)/2,1/2],[np.sqrt(2)/2,np.sqrt(2)/2],[1/2,np.sqrt(3)/2],[0,1],[-1/2,np.sqrt(3)/2],[-np.sqrt(2)/2,np.sqrt(2)/2],[-np.sqrt(3)/2,1/2],[-1,0],[-np.sqrt(3)/2,-1/2],[-np.sqrt(2)/2,-np.sqrt(2)/2],[-1/2,-np.sqrt(3)/2]]
    Z = 0
    p = np.zeros((len(qs),pfc.xRes,pfc.yRes))
    for q,n in zip(qs,range(len(qs))):
        eta = sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())
        # pfc.plot_field(abs(eta))
        # plt.show()
        p[n] = np.exp(abs(eta)/T)
        Z += p[n]
    
    q_net = np.zeros((2,pfc.xRes,pfc.yRes))
    for q,n in zip(qs,range(len(qs))):
        # pfc.plot_field(p[n]/Z)
        # plt.show()
        q_net[0] += q[0]*p[n]/Z
        q_net[1] += q[1]*p[n]/Z
    
    return q_net

q_net = orientation(pfc,0.01)

fig = plt.figure()
fig.subplots(1,2)

pfc.plot_field(pfc.psi,ax=fig.axes[0])

pfc.plot_vector_field(q_net,spacing=1,ax=fig.axes[1])
# pfc.plot_field(q_net[1])
plt.show()
            


# pfc.plot_field(pfc.psi)
# plt.show()

# fig = plt.figure()
# axs = fig.subplots(1,3)

# vlim = [0,0.02]
# pfc.plot_field(abs(demodulate(pfc,[1,0])),ax=axs[0],vlim=vlim)
# pfc.plot_field(abs(demodulate(pfc,[np.sqrt(2)/2,np.sqrt(2)/2])),ax=axs[1],vlim=vlim)
# pfc.plot_field(abs(demodulate(pfc,[0,1])),ax=axs[2],vlim=vlim)
# plt.show()


