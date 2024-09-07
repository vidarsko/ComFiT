import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import comfit as cf
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go



# bs = cf.BaseSystem(3,xRes=11,yRes=11,zRes=11)
# bs.plot_lib='plotly'
# vector_field = np.array([bs.z+bs.x*np.cos(bs.y/5), bs.z+bs.y*np.sin(bs.x/5)])
# fig = bs.plot_vector_field(vector_field)
# print(fig.layout['scene'])
# fig.show()

# plt.show()


# fig = plt.figure()

# #1D system
# bs = cf.BaseSystem(1,xRes=31)
# bs.plot_lib='plotly'

# 1D vector field
# vector_field = np.array([bs.x*np.cos(bs.x/5)])
# fig = bs.plot_vector_field(vector_field)
# fig.show()


# # 2D vector field needs work
# vector_field = np.array([bs.x*np.cos(bs.x/5), bs.x*np.sin(bs.x/5)])
# fig = bs.plot_vector_field(vector_field, spacing=2)
# fig.show()

# # # 3D vector field needs work
# vector_field = np.array([bs.x*np.cos(bs.x/5), bs.x*np.sin(bs.x/5), bs.x*np.cos(bs.x/5)])
# fig = bs.plot_vector_field(vector_field, spacing=3)
# fig.show()


# # 3D system
bs = cf.BaseSystem(3,xRes=11,yRes=11,zRes=11)
bs.plot_lib='plotly'

# # 1D vector field
# vector_field = np.array([bs.z+bs.x*np.cos(bs.y/5)])
# fig = bs.plot_vector_field(vector_field,spacing=3)
# fig.show()

# # 2D vector field
# vector_field = np.array([bs.z+ bs.x*np.cos(bs.y/5), bs.z + bs.y*np.sin(bs.x/5)])
# fig = bs.plot_vector_field(vector_field,spacing=5)
# fig.show()

# # 3D vector field
# vector_field = np.array([bs.z+ bs.x*np.cos(bs.y/5), bs.z + bs.y*np.sin(bs.x/5), -bs.z + bs.x*np.cos(bs.y/5)])
# fig = bs.plot_vector_field(vector_field,spacing=3)
# fig.show()

# #2D system
bs = cf.BaseSystem(2,xRes=31,yRes=31)
bs.plot_lib='plotly'

# # 1D vector field needs work
# vector_field = np.array([bs.x*np.cos(bs.y/5)])
# fig = bs.plot_vector_field(vector_field,spacing=3)
# fig.show()

# # 2D vector field
# vector_field = np.array([bs.x*np.cos(bs.y/5), bs.y*np.sin(bs.x/5)])
# fig = bs.plot_vector_field(vector_field,spacing=1)
# fig.show()

# # 3D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5), bs.y*np.sin(bs.x/5), bs.x*np.cos(bs.y/5)])
fig = bs.plot_vector_field(vector_field, spacing=3)
fig.show()


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


fig = plt.figure()

ax1 = fig.add_subplot(131)
bs = cf.BaseSystem(1,xRes=31)
field = bs.x**2
bs.plot_field(field,ax=ax1)

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