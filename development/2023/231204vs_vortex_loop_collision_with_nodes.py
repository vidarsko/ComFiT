import comfit as cf

import numpy as np

import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51,gamma=0)

#bec.conf_insert_vortex_ring(normal_vector=[-1,0,0],radius=10,position=[0.25*bec.xmax,0.5*bec.ymax,0.6*bec.zmax])

#bec.conf_insert_vortex_ring(normal_vector=[0,-1,0],radius=10)

bec.conf_insert_vortex_ring(normal_vector=[1,0,0],radius=10,
                            position=bec.rmid)
bec.evolve_relax(100)

# bec.plot_field(abs(bec.psi))

# plt.show()


N=10

for i in range(200):

    bec.evolve_dGPE(N-1)
    psi_old = bec.psi

    bec.evolve_dGPE(1)

    dt_psi =(bec.psi-psi_old)/(bec.dt)

    Dnodes = bec.calc_vortex_nodes(dt_psi)

    bec.plot_field(abs(bec.psi), colorbar=False,clims=[0,1],layer_values=[0.6])

    bec.plot_nodes(Dnodes)

    #theta = bec.calc_angle_field_vortex_ring(normal_vector=[0,0,1],radius=13)

    #bec.plot_angle_field(theta)

    cf.tool_zoom_plot_matplotlib(2)
    cf.tool_save_plot_matplotlib(i)

    #plt.draw()

    #plt.pause(0.03)

cf.tool_make_animation_movie(i)