import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(3,xRes=21,yRes=21,zRes=21)
bec.conf_insert_vortex_ring(normal_vector=[1,1,0])
bec.evolve_relax_BoseEinsteinCondensate(100)

bec.plot_field(abs(bec.psi))
plt.show()