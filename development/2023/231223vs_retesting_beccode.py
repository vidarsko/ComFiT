import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(2,xRes=51,yRes=51)
theta = bec.calc_angle_field_vortex_dipole()
bec.conf_insert_vortex_dipole(
    dipole_vector=[bec.xmax/3,0],
    dipole_position=bec.rmid)
bec.evolve_relax(100)
Dnodes = bec.calc_vortex_nodes()
bec.plot_field(abs(bec.psi),colorbar=False)
bec.plot_nodes(Dnodes)
plt.show()