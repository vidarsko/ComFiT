import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1,alpha=1)



nem.conf_insert_disclination_dipole(dipole_position=[32,12])

nem.conf_insert_disclination_dipole(dipole_position=[32,40])
nem.evolve_nematic_no_flow(10)



Q_prev = np.copy(nem.Q)
nem.evolve_nematic(10,"ETD4RK")

dt_Q = (nem.Q -Q_prev)/(10*nem.dt)

polarization = nem.calc_defect_polarization_field()
D = nem.calc_defect_density_nematic()
S,director =nem.calc_order_and_director()

Dnodes =nem.calc_disclination_nodes_nem(dt_Q =dt_Q)


ax =nem.plot_field_velocity_and_director(D,nem.u,director,colormap='spring')

nem.plot_disclination_nodes(Dnodes,ax=ax)


plt.show()
