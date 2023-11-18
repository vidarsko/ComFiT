import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=21,yRes=21,zRes=21)
bec.set_initial_condition_vortex_ring()
bec.evolve_relax_BEC(100)

bec.plot_field(abs(bec.psi))
plt.show()