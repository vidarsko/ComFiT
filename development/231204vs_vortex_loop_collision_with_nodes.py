import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=51,yRes=51,zRes=51,gamma=0)
bec.conf_insert_vortex_ring(normal_vector=[-1,0,0],radius=10,position=[0.25*bec.xmax,bec.ymid,bec.zmid])
#bec.conf_insert_vortex_ring(normal_vector=[0,-1,0],radius=10)
bec.conf_insert_vortex_ring(normal_vector=[1,0,0],radius=10,position=[0.75*bec.xmax,bec.ymid,bec.zmid])
bec.evolve_relax_BEC(100)
# bec.plot_field(abs(bec.psi))
# plt.show()

N=10
for i in range(400):
    psi_old = bec.psi
    bec.evolve_dGPE(N)
    dt_psi =(bec.psi-psi_old)/(N*bec.dt)
    Dnodes = bec.calc_vortex_nodes(dt_psi)
    bec.plot_field(abs(bec.psi), colorbar=False)
    bec.plot_vortex_nodes(Dnodes)
    #theta = bec.calc_angle_field_vortex_ring(normal_vector=[0,0,1],radius=13)
    #bec.plot_angle_field(theta)
    cf.tool_save_plot(i)
    #plt.draw()
    #plt.pause(0.03)
cf.tool_make_animation(i)