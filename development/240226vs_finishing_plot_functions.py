import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

ax1 = fig.add_subplot(231)
bs = cf.BaseSystem(1,xRes=31)
field = bs.x**2*np.exp(1j*bs.x/3)
bs.plot_complex_field(field,ax=ax1)

ax2 = fig.add_subplot(232)
bs = cf.BaseSystem(2,xRes=31,yRes=31)
field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
bs.plot_complex_field(field,ax=ax2,plot_method='phase_angle')

ax3 = fig.add_subplot(233, projection='3d')
bs = cf.BaseSystem(2,xRes=31,yRes=31)
field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
bs.plot_complex_field(field,ax=ax3,plot_method='3Dsurface')

ax5 = fig.add_subplot(235, projection='3d')
bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
bs.plot_complex_field(field,ax=ax5,plot_method='phase_angle')

ax6 = fig.add_subplot(236, projection='3d')
bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
bs.plot_complex_field(field,ax=ax6,plot_method='phase_blob')

plt.show()


# fig = plt.figure()

# ax1 = fig.add_subplot(131)
# bs = cf.BaseSystem(1,xRes=31)
# field = bs.x**2
# bs.plot_field(field,ax=ax1)

# ax2 = fig.add_subplot(132)
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# field = bs.x**2 + bs.y**2
# bs.plot_field(field,ax=ax2)

# ax3 = fig.add_subplot(133, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# field = bs.x**2 + bs.y**2 + bs.z**2
# bs.plot_field(field,ax=ax3)

# plt.show()


# bec.plot_field(np.abs(bec.psi))
# bec.plot_complex_field(bec.psi,axs[dim])




# bec = cf.BoseEinsteinCondensate(3,xRes=21,yRes=21,zRes=21)
# bec.conf_initial_condition_disordered()
# bec.evolve_relax(100)

# # bec.plot_field(np.abs(bec.psi))
# bec.plot_complex_field(bec.psi)

# plt.show()