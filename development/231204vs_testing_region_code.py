import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
#bec.conf_initial_condition_disordered()
# disk = bec.calc_region_disk([5,0],10)
#ball = bec.calc_region_ball([bec.xmid,0,0],10)
cylinder = bec.calc_region_cylinder([bec.xmid,bec.ymid,0],10,[0,1,1],10)
#psi_max_index = np.unravel_index(np.argmax(np.abs(bec.psi)), bec.psi.shape)
#print(psi_max_index[1])
# bec.plot_field(disk)
bec.plot_field(cylinder)
plt.show()