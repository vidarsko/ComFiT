import comfit as cf

import comfit as cf
import matplotlib.pyplot as plt


bec = cf.BEC(2,xRes=201,yRes=201)
bec.set_initial_condition_disordered()
bec.evolve_relax_BEC(200)

# rho = bec.calc_vortex_density_singular()
#
# bec.plot_field(rho)
# plt.show()

nodes = bec.calc_vortex_nodes()
bec.plot_vortices(nodes)
plt.show()
