import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt



pfc = cf.PhaseFieldCrystal2DTriangular(30,20)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)

for n in range(100):
    pfc.evolve_PFC(100)
    pfc.plot_field(pfc.psi, grid=False, colorbar=False)
    cf.tool_save_plot(n)

cf.tool_make_animation_gif(n)



# qm = cf.QuantumMechanics(1,xlim=[-10,100],dx=0.1)

# qm.V_ext = qm.calc_Gaussian(position=30,width=10,top=1)

# qm.conf_initial_condition_Gaussian(position=0,width=1,initial_velocity=1)

# ymax = np.max(np.abs(qm.psi)**2)

# for n in range(100):
#     qm.plot(ylim=ymax)
#     cf.tool_save_plot(n)
#     qm.evolve_schrodinger(5)

# cf.tool_make_animation_gif(n)



# qm.plot()
# plt.show()

# a=5
# b=10
# interval = qm.calc_region_interval(a,b)
# prob = qm.calc_integrate_field(np.abs(qm.psi)**2,interval)

# print("The probability of measiring the particle in the interval (",a,",",b,") is",prob)

# qm.evolve_schrodinger(45)
# qm.plot()
# plt.show()

# a=5
# b=10
# interval = qm.calc_region_interval(a,b)
# prob = qm.calc_integrate_field(np.abs(qm.psi)**2,interval)

# print("The probability of measiring the particle in the interval (",a,",",b,") is",prob)

# import mayavi.mlab as mlab

# bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
# f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
# bs.plot_field(f,plotting_lib='mayavi', number_of_layers=3)
# mlab.show()


# import mayavi.mlab as mlab

# bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
# f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
# bs.plot_field_in_plane(f, position=[0,0,0], normal_vector=[1,1,1],plotting_lib='mayavi')
# mlab.show()





# import matplotlib.pyplot as plt



# bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
# f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
# bs.plot_field_in_plane(f, position=[0,0,0], normal_vector=[1,1,0])
# plt.show()

# bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
# f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
# bs.plot_field(f,number_of_layers=3)
# plt.show()

# bs = cf.BaseSystem(2, xlim=[-3,3], ylim=[-3,3])
# f = np.sin(bs.x+bs.y)
# bs.plot_field(f)
# plt.show()