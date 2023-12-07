import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101)
bec.conf_insert_vortex_dipole()
bec.evolve_relax_BoseEinsteinCondensate(300)
#Dnodes = bec.calc_vortex_nodes()

psi_old = bec.psi
N=10
bec.evolve_dGPE(N)
dt_psi = (bec.psi-psi_old)/(N*bec.dt)

# psi = [np.real(bec.psi),np.imag(bec.psi)]
# dt_psi = [np.real(dt_psi),np.imag(dt_psi)]
Dnodes = bec.calc_vortex_nodes(dt_psi)
print(Dnodes)

#defect_field = bec.calc_defect_density(psi)
#velocity_field = bec.calc_defect_velocity_field(psi,dt_psi)

#bec.plot_vector_field(velocity_field, step = 2)

#bec.plot_field(velocity_field[1])

bec.plot_vortex_nodes(Dnodes)
plt.show()