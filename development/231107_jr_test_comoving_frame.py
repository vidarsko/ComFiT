import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BEC(2,xRes=512,yRes=256,dx=0.5,dy=0.5,gamma=0,dt=0.01)

bec.V_ext = bec.gaussian_stirring_potential(4,10,[bec.xmid,bec.ymid])
bec.set_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(100)

bec.set_spatialy_varying_gamma(wx=100,wy=50)

plt.imshow(bec.V_ext)
plt.colorbar()
plt.show()

for i in range(10):
    bec.evolve_comoving_dGPE(1000,0.35)
    plt.imshow(np.abs(bec.psi)**2)
    plt.show()