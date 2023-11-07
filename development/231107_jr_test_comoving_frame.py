import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BEC(2,xRes=256,yRes=128,gamma=0,dt=0.01)

bec.V_ext = bec.gaussian_stirring_potential(4,10,[bec.xmid+50,bec.ymid])

bec.set_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(1000)

#bec.set_spatialy_varying_gamma(wx=100,wy=50)

k2 = bec.calc_k2()


#plt.imshow(bec.gamma)
#plt.colorbar()
#plt.show()

for i in range(10):
    bec.evolve_comoving_dGPE(1000,0.5)
   # bec.evolve_relax_BEC(10000)
    plt.imshow(np.abs(bec.psi)**2)
    plt.colorbar()
    plt.show()