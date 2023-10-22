import comfit as cf

import comfit as cf
import matplotlib.pyplot as plt


bec = cf.BEC(2,xRes=101,yRes=101)
bec.set_initial_condition_disordered()
bec.evolve_relax_BEC(100)

# rho = bec.calc_vortex_density_singular()
#
# bec.plot_field(rho)
# plt.show()

nodes = bec.calc_vortex_nodes()
print(nodes)

