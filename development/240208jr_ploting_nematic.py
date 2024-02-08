import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(3,xRes=31,yRes=31,zRes = 31,dx=1,dy=1,dt=0.1,alpha=-.5,C=1)



nem.conf_initial_condition_disordered(noise_strength=1.0)
nem.evolve_nematic_no_flow(20,method="ETD4RK")


nem.evolve_nematic(600,method="ETD4RK")

S,n = nem.calc_order_and_director()
nem.plot_field(S, plotting_lib='mayavi', number_of_layers=2)
mlab.show()