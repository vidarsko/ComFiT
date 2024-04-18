import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

pfc = cf.PhaseFieldCrystal2DTriangular(20,13)
eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole())
# q0=[np.sqrt(2)/2,np.sqrt(2)/2]
# q0 = [0,1]sd
# q0 = [np.sqrt(3)/2,1/2]
# pfc.psi = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))

pfc.psi = np.random.rand(pfc.xRes,pfc.yRes)-0.5

pfc.psi_f = np.fft.fft2(pfc.psi)

pfc.a0 = 2*np.pi

pfc.evolve_PFC(500)

def demodulate(pfc,q):
    return sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())


def orientation(pfc,T):
    N=8
    qs = [[np.cos(theta),np.sin(theta)] for theta in np.linspace(0,np.pi-(np.pi)/N,N)]
    Z = 0
    p = np.zeros((len(qs),pfc.xRes,pfc.yRes))
    for q,n in zip(qs,range(len(qs))):
        eta = sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())
        # pfc.plot_field(abs(eta))
        # plt.show()
        p[n] = np.exp(abs(eta)/T)
        Z += p[n]
    
    Q_net = np.zeros((2,pfc.xRes,pfc.yRes))
    for q,n in zip(qs,range(len(qs))):
        # print(q[0]*q[0]-1/2)
        # print(q[0]*q[1])
        # pfc.plot_field(p[n]/Z)
        # plt.show()
        prob = p[n]/Z
        Q_net[0] += prob*(q[0]*q[0]-1/2)
        Q_net[1] += prob*q[0]*q[1]
    
    return Q_net

Q_net = orientation(pfc,0.1)
theta = np.arctan2(Q_net[1],Q_net[0])
q = [np.cos(theta/2),np.sin(theta/2)]


fig = plt.figure()
fig.subplots(1,1)

pfc.plot_field(pfc.psi,ax=fig.axes[0])

pfc.plot_vector_field(q,spacing=1,ax=fig.axes[0])
# pfc.plot_field(Q_net[0],ax=fig.axes[2])
# pfc.plot_field(Q_net[1],ax=fig.axes[3])
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


