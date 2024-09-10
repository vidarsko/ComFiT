import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1)



nem.conf_initial_condition_disordered(noise_strength=0.2)

nem.conf_active_channel(20,7)

nem.evolve_nematic(4000,"ETD4RK")

D = nem.calc_defect_density_nematic()
director =nem.calc_director()
ax= nem.plot_field_velocity_and_director(D,nem.u,director)
plt.show()