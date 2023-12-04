import comfit as cf
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=31,yRes=31,zRes=31)
bec.conf_insert_vortex_ring()
bec.evolve_relax_BEC(300)

Dnodes = bec.calc_vortex_nodes()
print(Dnodes)
bec.plot_vortex_nodes(Dnodes)
plt.show()

