import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem =cf.NematicLiquidCrystal(3,xRes=32,yRes=32,zRes=32,dx=1,dt=0.1,alpha=-0.8,C=1)



nem.conf_initial_disclination_lines()

nem.evolve_nematic_no_flow( 170,method="ETD2RK")
Q_prev = np.copy(nem.Q)
nem.evolve_nematic_no_flow(1)

dt_Q = (nem.Q - Q_prev)/(nem.dt)

nodes = nem.calc_disclination_nodes_nem(dt_Q)

#print(nodes[3]['velocity'])

fig,axes =nem.plot_nodes(nodes)


fig.show()
