import comfit as cf
import numpy as np
import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=61,yRes=61,zRes=61)
bec.conf_insert_vortex_ring()
#bec.evolve_relax(100)

bec.plot_field_in_plane(np.angle(bec.psi))
plt.show()

