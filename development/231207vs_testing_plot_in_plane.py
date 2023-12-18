import comfit as cf
import numpy as np
import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31,
                                dx=2,dy=2,dz=2)
bec.conf_insert_vortex_ring()
bec.evolve_relax(100)
#bec.conf_initial_condition_disordered()

bec.plot_field_in_plane(np.abs(bec.psi),normal_vector=[1,1,1],
                         position=[bec.xmid+0.01,bec.ymid+0.03,bec.zmid+0.32])
# bec.plot_field_in_plane(np.abs(bec.psi),normal_vector=[1,1,1],
#                         position=[bec.xmid,bec.ymid,bec.zmid])
plt.show()

