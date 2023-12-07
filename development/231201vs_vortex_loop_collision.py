import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51,gamma=0)
bec.conf_insert_vortex_ring(normal_vector=[0,-1,0],radius=10,position=[bec.xmid,bec.ymax,bec.zmax/2])
#bec.conf_insert_vortex_ring(normal_vector=[0,-1,0],radius=10)
bec.conf_insert_vortex_ring(normal_vector=[-1,0,0],radius=10,position=[bec.xmax,bec.ymid,bec.zmax/2])
bec.evolve_relax_BoseEinsteinCondensate(100)
# bec.plot_field(abs(bec.psi))
# plt.show()

for i in range(20):
    bec.evolve_dGPE(10)
    bec.plot_field(abs(bec.psi))
    #theta = bec.calc_angle_field_vortex_ring(normal_vector=[0,0,1],radius=13)
    #bec.plot_angle_field(theta)
    cf.tool_save_plot(i)
    #plt.draw()
    #plt.pause(0.03)
cf.tool_make_animation(i)