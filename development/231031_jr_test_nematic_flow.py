import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.Nematic(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1)



nem.conf_initial_condition_disordered(noise_strength=4)
nem.evolve_nematic_no_flow(10)

nem.evolve_nematic(4000)
D = nem.calc_defect_density_nematic()
director =nem.calc_director()
ax= nem.plot_field_velocity_and_director(D,nem.u,director)

plt.show()