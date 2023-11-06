import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.nematic(2,xRes=32,yRes=32,dx=1,dy=1,dt=0.01)

nem.set_initial_condition_disordered(noise_strength=0.5)
nem.evolve_nematic_no_flow(5)


nem.evolve_nematic(10)
D = nem.calc_defect_density_nematic()
director =nem.calc_director()
ax= nem.plot_field_velocity_and_director(D,nem.u,director)

plt.show()