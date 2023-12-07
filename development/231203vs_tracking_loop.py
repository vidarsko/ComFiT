import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51)
bec.conf_insert_vortex_ring()
bec.evolve_relax_BoseEinsteinCondensate(300)
bec.plot_field(np.abs(bec.psi)**2)
