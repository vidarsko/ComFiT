import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201,gamma=0)
bec.conf_initial_condition_disordered()
bec.evolve_relax_BoseEinsteinCondensate(300)

for n in range(300):
    bec.evolve_dGPE(20)
    bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False,
                   clims=[0,1])
    cf.tool_save_plot(n)

cf.tool_make_animation(n)


# bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201)
# bec.conf_initial_condition_disordered()
#
# for n in range(300):
#     bec.evolve_relax_BoseEinsteinCondensate(1)
#     bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False,
#                    clims=[0,1])
#     cf.tool_save_plot(n)
#
# cf.tool_make_animation(n)


bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201,gamma=0)
bec.conf_initial_condition_disordered()
bec.evolve_relax_BoseEinsteinCondensate(300)

for n in range(300):
    bec.evolve_dGPE(20)
    bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False,
                   clims=[0,1])
    cf.tool_save_plot(n)

cf.tool_make_animation(n)


# bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201,gamma=0)
# bec.conf_insert_vortex_dipole(dipole_vector=[0.5*bec.xmax,0])
# bec.conf_insert_vortex_dipole(dipole_position=[0.625*bec.xmax,bec.ymid],dipole_vector=[0.25*bec.xmax,0])
# bec.evolve_relax_BoseEinsteinCondensate(0)
# #bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False, clims=[0,1])
# bec.plot_angle_field(np.angle(bec.psi))
# #cf.tool_save_plot(0)
# plt.show()


# bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51,gamma=0)
# bec.conf_initial_condition_disordered()
# for n in range(300):
#     bec.evolve_relax_BoseEinsteinCondensate(1)
#     bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False,
#                    clims=[0,1])
#     cf.tool_save_plot(n)
#
# cf.tool_make_animation(n)

# bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51,gamma=0)
# bec.conf_initial_condition_disordered()
# bec.evolve_relax_BoseEinsteinCondensate(300)
#
# for n in range(300):
#     bec.evolve_dGPE(20)
#     bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False,
#                    clims=[0,1],layer_values=[0.6])
#     cf.tool_save_plot(n)
#
# cf.tool_make_animation(n)


# bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201,gamma=0)
# bec.conf_insert_vortex_dipole()
# # bec.conf_insert_vortex_dipole(dipole_position=[0.625*bec.xmax,bec.ymid],dipole_vector=[0.25*bec.xmax,0])
# bec.evolve_relax_BoseEinsteinCondensate(100)
# #bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False, clims=[0,1])
# # bec.plot_angle_field(np.angle(bec.psi))
# rho = bec.calc_vortex_density()
# bec.plot_field(rho)
# #cf.tool_save_plot(0)
# plt.show()


bec = cf.BoseEinsteinCondensate(2,xRes=201,yRes=201)
bec.conf_initial_condition_disordered()
bec.evolve_relax_BoseEinsteinCondensate(300)
rho = bec.calc_vortex_density()
bec.plot_field(rho)
cf.tool_save_plot(0)
bec.plot_angle_field(np.angle(bec.psi))
cf.tool_save_plot(1)
bec.plot_field(abs(bec.psi),colormap='winter',cmap_symmetric=False, clims=[0,1])
# plt.show()y
