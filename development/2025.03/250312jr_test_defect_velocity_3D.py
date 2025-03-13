import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem =cf.NematicLiquidCrystal(3,xRes=32,yRes=32,zRes=32,dx=1,dt=0.1,alpha=-0.8)



nem.conf_initial_disclination_lines()

nem.evolve_nematic_no_flow(90,method="ETD4RK")
Q_prev = np.copy(nem.Q)
nem.evolve_nematic_no_flow(5)

dt_Q = (nem.Q - Q_prev)/(5*nem.dt)

nodes = nem.calc_disclination_nodes_nem(dt_Q=dt_Q)

#print(nodes[3]['velocity'])

fig,axes =nem.plot_nodes(nodes)

fig.show()
