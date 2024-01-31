import comfit as cf

import numpy as np

import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=101,yRes=101,zRes=101,gamma=0)

# print(bec.psi)

# bec.conf_insert_vortex_ring(normal_vector=[0,1,0],radius=20)

# bec.conf_insert_vortex_ring(normal_vector=[0,-1,0],radius=10)
bec.conf_initial_condition_disordered()


bec.evolve_relax(1000)


for i in range(1000):

    bec.evolve_dGPE(10)
    bec.plot_field(abs(bec.psi))

    #theta = bec.calc_angle_field_vortex_ring(normal_vector=[0,0,1],radius=13)

    #bec.plot_angle_field(theta)
    cf.tool_save_plot(i)

    #plt.draw()

    #plt.pause(0.03)

cf.tool_make_animation_movie(i)