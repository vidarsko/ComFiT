import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=31,yRes=31,zRes=31)
bec.conf_insert_vortex_ring()
bec.evolve_relax_BEC(300)


rho = bec.calc_defect_density([np.real(bec.psi),np.imag(bec.psi)])

rho_abs = np.sqrt(rho[0]**2+rho[1]**2+rho[2]**2)

#bec.plot_field(rho_abs)
bec.plot_vector_field(rho)
cf.tool_zoom_plot(2)
#cf.tool_export_rotating_plot()
plt.show()