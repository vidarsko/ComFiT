import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

N=16
T=0.1

pfc = cf.PhaseFieldCrystal2DTriangular(40,26)
eta = np.exp(1j*pfc.calc_angle_field_vortex_dipole())
q0 = [1,0]
# q0=[np.sqrt(2)/2,np.sqrt(2)/2]
# q0 = [0,1]
# q0 = [np.sqrt(3)/2,1/2]
theta=np.pi/N/2
q0 = [np.cos(theta),np.sin(theta)]
pfc.psi = np.real(eta*np.exp(1j*(q0[0]*pfc.x + q0[1]*pfc.y)))

# pfc.psi = np.random.rand(pfc.xRes,pfc.yRes)-0.5

pfc.psi_f = np.fft.fft2(pfc.psi)

pfc.a0 = 2*np.pi

pfc.evolve_PFC(500)

def demodulate(pfc,q):
    return sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())


def orientation(pfc,T):
    qs = [[np.cos(theta),np.sin(theta)] for theta in np.linspace(0,np.pi-(np.pi)/N,N)]
    Z = 0
    p = np.zeros((len(qs),pfc.xRes,pfc.yRes))
    eta = np.zeros((len(qs),pfc.xRes,pfc.yRes),dtype=complex)
    for q,n in zip(qs,range(len(qs))):
        eta[n] = sp.fft.ifftn(sp.fft.fftn(pfc.psi*np.exp(-1j*q[0]*pfc.x - 1j*q[1]*pfc.y))*pfc.calc_Gaussian_filter_f())
        # pfc.plot_field(abs(eta))
        # plt.show()
        p[n] = np.exp(abs(eta[n])/T)
        Z += p[n]
    
    Q_net = np.zeros((2,pfc.xRes,pfc.yRes))
    eta_net = np.zeros((pfc.xRes,pfc.yRes),dtype=complex)
    eta_mag = np.zeros((pfc.xRes,pfc.yRes))
    eta_phase = np.ones((pfc.xRes,pfc.yRes),dtype=complex)
    for q,n in zip(qs,range(len(qs))):
        # print(q[0]*q[0]-1/2)
        # print(q[0]*q[1])
        # pfc.plot_field(p[n]/Z)
        # plt.show()
        # pfc.plot_complex_field(eta[n])
        # plt.show()
        prob = p[n]/Z
        Q_net[0] += prob*(q[0]*q[0]-1/2)
        Q_net[1] += prob*q[0]*q[1]
        eta_net += prob*eta[n]
        eta_mag += prob*abs(eta[n])
        eta_phase *= np.exp(1j*prob*np.angle(eta[n]))
    eta_net2 = eta_phase*eta_mag
    
    return Q_net, eta_net, eta_net2

Q_net, eta_net,eta_net2 = orientation(pfc,T)

def coarse_grain(field,width):
    return sp.fft.ifftn(sp.fft.fftn(field)*pfc.calc_Gaussian_filter_f(a0=width))

Q_net[0] = coarse_grain(Q_net[0],5*pfc.a0)
Q_net[1] = coarse_grain(Q_net[1],5*pfc.a0)
theta = np.arctan2(Q_net[1],Q_net[0])

q = [np.cos(theta/2),np.sin(theta/2)]

# psi_new = np.real(eta*np.exp(1j*(q[0]*pfc.x + q[1]*pfc.y)))
# # pfc.plot_complex_field(psi_new)
# pfc.plot_field(psi_new)
# plt.show()

eta_supreme = demodulate(pfc,q)
# eta_supreme = demodulate(pfc,[1,0])
# angle=0.1
# eta_supreme = demodulate(pfc,[np.cos(angle),np.sin(angle)])

# pfc.plot_vector_field(q,spacing=1)
# plt.show()

fig = plt.figure()
fig.subplots(2,2)

pfc.plot_field(pfc.psi,ax=fig.axes[0])

# pfc.plot_complex_field(eta_net,ax=fig.axes[1])
# fig.axes[1].text(0.1, 0.9, f'max eta-value: {np.max(abs(eta_net)):.2f}')

# pfc.plot_complex_field(eta_net2,ax=fig.axes[2])
# fig.axes[2].text(0.1, 0.9, f'max eta-value: {np.max(abs(eta_net2)):.2f}')

pfc.plot_vector_field(q,spacing=15,ax=fig.axes[1])
pfc.plot_angle_field(np.arctan2(q[1],q[0]),ax=fig.axes[2])
# pfc.plot_complex_field(eta_supreme,ax=fig.axes[3])
# fig.axes[3].text(0.1, 0.9, f'max eta-value: {np.max(abs(eta_supreme)):.2f}')
pfc.plot_field(abs(eta_supreme),ax=fig.axes[3])
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


