import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem =cf.NematicLiquidCrystal(3,xRes=64,yRes=64,zRes=64,dt=0.1,alpha=-1.2,C=1)

#nem.conf_initial_condition_ordered(noise_strength=0.1)

nem.conf_initial_disclination_lines()


nem.evolve_nematic_no_flow(100,method="ETD4RK")
#nem.evolve_nematic(2000,method="ETD4RK")

S,n = nem.calc_order_and_director()
#D = nem.calc_disclination_density_nematic()

omega,Omega,T,trD = nem.calc_disclination_density_decoupled()
#fig,ax =nem.plot_field(omega,number_of_layers=2,vlim_symmetric=False)
#nem.plot_vector_field(Omega)

nodes = nem.calc_disclination_nodes_nem()
#print(nodes)
nem.plot_nodes(nodes)

plt.show()