import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.nematic(2,xRes=101,yRes=101,dt=0.1)

nem.set_initial_condition_disordered(noise_strength=.5)
nem.evolve_nematic_no_flow(200)

S =nem.calc_S()

plt.imshow(S)
plt.colorbar()
plt.show()