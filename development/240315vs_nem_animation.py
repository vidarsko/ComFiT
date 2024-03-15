# Task 1 + 2 for animation see e.g 2.2 - BoseEinsteinCondensate tutorial: Time dependent potentials
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1,alpha=1)



nem.conf_insert_disclination_dipole()
nem.evolve_nematic_no_flow(10)



Q_prev = np.copy(nem.Q)

for n in range(200):
    
    nem.evolve_nematic(3,"ETD4RK")
    dt_Q = (nem.Q -Q_prev)/(10*nem.dt)

    polarization = nem.calc_defect_polarization_field()
    D = nem.calc_defect_density_nematic()
    S,director =nem.calc_order_and_director()

    Dnodes =nem.calc_disclination_nodes_nem(dt_Q =dt_Q)

    plt.clf()
    ax =nem.plot_field_velocity_and_director(D,nem.u,director,colormap='spring',cmin=-0.15, cmax=0.15)

    nem.plot_disclination_nodes(Dnodes,ax=ax)
    cf.tool_save_plot(n)

cf.tool_make_animation_gif(n)