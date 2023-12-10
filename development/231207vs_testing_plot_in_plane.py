import comfit as cf
import numpy as np
import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=41,yRes=41,zRes=41)
bec.conf_insert_vortex_ring()
bec.evolve_relax(100)

bec.plot_field_in_plane(np.abs(bec.psi),normal_vector=[1,1,1],
                        position=[bec.xmid,bec.ymid,bec.zmid+5])
plt.show()

