import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1)

nem.conf_insert_disclination_dipole()
nem.evolve_nematic_no_flow(10,method = 'ETD2RK')

S,n = nem.calc_order_and_director()

nem.evolve_nematic_no_flow(10)

Q_prev = np.copy(nem.Q)
nem.evolve_nematic(10,"ETD4RK")

dt_Q = (nem.Q -Q_prev)/(10*nem.dt)

polarization = nem.calc_disclination_polarization_field()
D = nem.calc_disclination_density_nematic()
S,director =nem.calc_order_and_director()

Dnodes =nem.calc_disclination_nodes_nem(dt_Q=dt_Q)

fig = nem.plot_field_velocity_and_director(D,nem.u,director)
# nem.plot_disclination_nodes(Dnodes,fig=fig)
fig.show()


# nem.plot_lib='matplotlib'
# nem.plot_field_velocity_and_director(D,nem.u,director)
# plt.show()

# 