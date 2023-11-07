import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BEC(3,xRes=64,yRes=64,zRes=64,gamma=0,dt=0.1)

bec.V_ext = bec.gaussian_stirring_potential(4,0.8,[bec.xmid,bec.ymid,bec.zmid])

bec.set_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(100)

bec.set_spatialy_varying_gamma(wx=20,wy=20,wz=20)

k2 = bec.calc_k2()


#plt.imshow(bec.gamma)
#plt.colorbar()
#plt.show()
bec.psi += (0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)+ 0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)*(1j))*np.abs(bec.psi)**2

for i in range(10):
    bec.evolve_comoving_dGPE(1000,0.4)
   # bec.evolve_relax_BEC(10000)
    bec.plot_field(np.abs(bec.psi)**2)
    plt.show()