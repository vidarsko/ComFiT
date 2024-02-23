import numpy as np
import matplotlib.pyplot as plt
import comfit as cf
from mayavi import mlab

nem = cf.NematicLiquidCrystal(3,xRes=31,yRes=31,zRes = 31,dx=1,dy=1,dt=0.1,alpha=-1.5,C=1)



nem.conf_initial_condition_disordered(noise_strength=1)

#nem.conf_initial_disclination_line()

nem.evolve_nematic_no_flow(5,method="ETD4RK")



#for i in range(600):

#    nem.evolve_nematic(10)
#    S,n = nem.calc_order_and_director()
#    nem.plot_field(S,number_of_layers=2,vlim_symmetric=False)
#    plt.draw()
#    plt.pause(0.01)
#    cf.tool_save_plot(i)
#    plt.clf()
cf.tool_make_animation_gif(539)

plt.show()