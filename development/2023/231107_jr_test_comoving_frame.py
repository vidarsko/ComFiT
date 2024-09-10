import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BoseEinsteinCondensate(3,xRes=64,yRes=64,zRes=64,gamma=0,dt=0.1)

bec.V0 = bec.calc_gaussian_stirring_potential(2,4,[bec.xmid,bec.ymid,bec.zmid])

bec.conf_initial_condition_Thomas_Fermi()
bec.evolve_relax(100)

bec.conf_dissipative_frame(wx=26,wy=26,wz=26)



#print(np.min(k2))
#plt.imshow(bec.gamma)
#plt.colorbar()
#plt.show()
bec.psi += (0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)+ 0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)*(1j))*np.abs(bec.psi)**2

for i in range(1):
    bec.evolve_comoving_dGPE(1000,0.8)
   # bec.evolve_relax(10000)
    bec.plot_field(np.abs(bec.psi)**2)
    plt.show()