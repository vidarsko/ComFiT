import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

bec = cf.BoseEinsteinCondensate(3,xRes=64,yRes=64,zRes=64,gamma=0,dt=0.05)

pot = bec.calc_gaussian_stirring_potential(3.5,4,[bec.xmid,bec.ymid,bec.zmid])

bec.conf_external_potential(pot, additive=False)

bec.conf_initial_condition_Thomas_Fermi()
bec.evolve_relax(100)

bec.conf_dissipative_frame(wx=25,wy=25,wz=25)



## add noise to break symmetry
bec.psi += (0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)+ 0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)*(1j))*np.abs(bec.psi)**2
bec.psi_f = np.fft.fftn(bec.psi)
vel_x = 0.8

N = 600



for n in range(N):
    bec.evolve_comoving_dGPE(30,vel_x,method='ETD4RK')
      

    bec.plot_complex_field(bec.psi)
  
    cf.tool_save_plot(n)


cf.tool_make_animation_gif(n)