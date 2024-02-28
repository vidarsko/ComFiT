import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()


bs = cf.BaseSystem(1,xRes=31)

# 1D vector field
ax1 = fig.add_subplot(331)
vector_field = np.array([bs.x**2*np.cos(bs.x/5)])
bs.plot_vector_field(vector_field,ax=ax1)

# 2D vector field
ax2 = fig.add_subplot(332, projection='3d')
vector_field = np.array([bs.x**2*np.cos(bs.x/5), bs.x**2*np.sin(bs.x/5)])
bs.plot_vector_field(vector_field,ax=ax2, spacing=3)

# 3D vector field
ax3 = fig.add_subplot(333, projection='3d')
vector_field = np.array([bs.x**2*np.cos(bs.x/5), bs.x**2*np.sin(bs.x/5), bs.x**2*np.cos(bs.x/5)])
bs.plot_vector_field(vector_field,ax=ax3, spacing=3)

plt.show()


# import comfit as cf
# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure()

# ax1 = fig.add_subplot(121, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# angle_field = np.mod((bs.x + 2*bs.y + 3*bs.z)/5,2*np.pi)-np.pi
# bs.plot_angle_field_in_plane(angle_field, ax=ax1)

# ax2 = fig.add_subplot(122, projection='3d')
# bs.plot_angle_field_in_plane(angle_field, ax=ax2, normal_vector=[0,0,1],position=[10,10,10])

# plt.show()

# import comfit as cf
# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure()

# ax1 = fig.add_subplot(131)
# bs = cf.BaseSystem(1,xRes=31)
# angle_field = np.mod((bs.x)/5,2*np.pi)-np.pi
# bs.plot_angle_field(angle_field,ax=ax1)

# ax2 = fig.add_subplot(132)
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# angle_field = np.mod((bs.x + 2*bs.y)/5,2*np.pi)-np.pi
# bs.plot_angle_field(angle_field,ax=ax2)

# ax3 = fig.add_subplot(133, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# angle_field = np.mod((bs.x + 2*bs.y + 3*bs.z)/5,2*np.pi)-np.pi
# bs.plot_angle_field(angle_field,ax=ax3)

# plt.show()


# fig = plt.figure()

# ax1 = fig.add_subplot(231)
# bs = cf.BaseSystem(1,xRes=31)
# field = bs.x**2*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax1)

# ax2 = fig.add_subplot(232)
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax2,plot_method='phase_angle')

# ax3 = fig.add_subplot(233, projection='3d')
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax3,plot_method='3Dsurface')

# ax5 = fig.add_subplot(235, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax5,plot_method='phase_angle')

# ax6 = fig.add_subplot(236, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax6,plot_method='phase_blob')

# plt.show()

# fig = plt.figure()

# ax1 = fig.add_subplot(121, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# complex_field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.y/3)
# bs.plot_complex_field_in_plane(complex_field, ax=ax1)

# ax2 = fig.add_subplot(122, projection='3d')
# bs.plot_complex_field_in_plane(complex_field, ax=ax2, normal_vector=[0,0,1],position=[10,10,10])

# plt.show()

# ax1 = fig.add_subplot(131
# bs = cf.BaseSystem(1,xRes=31)
# vector_field = np.array([bs.x**2, bs.x**3])
# bs.plot_vector_field(vector_field,ax=ax1)
# plt.show()


# fig = plt.figure()

# ax1 = fig.add_subplot(231)
# bs = cf.BaseSystem(1,xRes=31)
# field = bs.x**2*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax1)

# ax2 = fig.add_subplot(232)
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax2,plot_method='phase_angle')

# ax3 = fig.add_subplot(233, projection='3d')
# bs = cf.BaseSystem(2,xRes=31,yRes=31)
# field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax3,plot_method='3Dsurface')

# ax5 = fig.add_subplot(235, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax5,plot_method='phase_angle')

# ax6 = fig.add_subplot(236, projection='3d')
# bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
# field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
# bs.plot_complex_field(field,ax=ax6,plot_method='phase_blob')

# plt.show()


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