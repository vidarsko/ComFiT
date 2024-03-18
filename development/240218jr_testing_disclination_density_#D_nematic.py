import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem =cf.NematicLiquidCrystal(3,xRes=32,yRes=32,zRes=32,dt=0.1,alpha=-1.5,C=1)

nem.conf_initial_condition_ordered(noise_strength=0.1)

#nem.conf_initial_disclination_line()


#nem.evolve_nematic_no_flow(5,method="ETD4RK")
nem.evolve_nematic(2500,method="ETD4RK")

S,n = nem.calc_order_and_director()
D = nem.calc_disclination_density_nematic()
omega = nem.calc_disclination_density_decoupling()
nem.plot_field(omega,number_of_layers=2,vlim_symmetric=False)

plt.show()