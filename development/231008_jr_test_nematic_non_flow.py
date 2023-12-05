import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.Nematic(2,xRes=101,yRes=101,dx=1,dy=1,dt=0.1)


nem.conf_initial_condition_disordered(noise_strength=4)

nem.evolve_nematic_no_flow(50)


S =nem.calc_S()

plt.imshow(S)
plt.colorbar()
plt.show()