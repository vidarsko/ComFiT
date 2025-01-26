import numpy as np
import matplotlib.pyplot as plt
import comfit as cf
from mayavi import mlab

nem = cf.NematicLiquidCrystal(2,xRes=65,yRes=65,dt=0.1)


#nem.conf_initial_condition_disordered(noise_strength=1)

nem.conf_insert_disclination_dipole()

nem.evolve_nematic_no_flow(5,method="ETD4RK")

nem.evolve_nematic(20)

s,n = nem.calc_order_and_director()

nem.plot_field_velocity_and_director(s,nem.u,n)

plt.show()
