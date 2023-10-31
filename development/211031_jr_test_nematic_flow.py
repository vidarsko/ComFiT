import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.nematic(2,xRes=101,yRes=101,dx=1,dy=1,dt=0.1)

nem.set_initial_condition_disordered(noise_strength=4)
nem.evolve_nematic_no_flow(50)


nem.evolve_nematic(20000)
D = nem.calc_defect_density_nematic()
director =nem.calc_director()
ax= nem.plot_field_velocity_and_director(D,nem.u,director)

plt.show()