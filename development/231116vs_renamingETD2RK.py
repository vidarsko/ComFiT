import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=51,yRes=51)
bec.set_initial_condition_disordered()
bec.evolve_relax_BEC(100)
rho = bec.calc_vortex_density()

bec.plot_field(rho)
plt.show()