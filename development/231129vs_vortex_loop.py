import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=31,yRes=31,zRes=31,gamma=0)
bec.set_initial_condition_vortex_ring(normal_vector=[1,0,0],radius=11)
bec.evolve_relax_BEC(100)
ax=None
for i in range(1000):
    bec.evolve_dGPE(10)
    ax = bec.plot_field(abs(bec.psi),ax=ax,colorbar=False)
    #theta = bec.calc_angle_field_vortex_ring(normal_vector=[0,0,1],radius=13)
    #bec.plot_angle_field(theta)
    cf.tool_save_plot(i)

cf.tool_make_animation(i)