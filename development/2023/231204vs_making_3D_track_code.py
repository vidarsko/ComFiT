import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
bec.conf_insert_vortex_ring()
bec.evolve_relax(300)

psi_vec = [np.real(bec.psi),np.imag(bec.psi)]
psi_old = bec.psi
bec.evolve_dGPE(1)
dt_psi =(bec.psi-psi_old)/bec.dt
dt_psi_vec = [np.real(dt_psi),np.imag(dt_psi)]

# velocity_field = bec.calc_defect_velocity_field(psi_vec,dt_psi_vec)
# bec.plot_field(velocity_field[2],number_of_layers=2)
Dnodes = bec.calc_vortex_nodes(dt_psi)
# print(Dnodes)
bec.plot_field(abs(bec.psi),colorbar=False)
bec.plot_nodes(Dnodes)
plt.show()

