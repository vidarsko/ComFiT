import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem =cf.NematicLiquidCrystal(2,xRes=101,yRes=101,dx=1,dt=0.1,alpha=-0.8)

nem.conf_initial_condition_ordered(noise_strength=0.1)

#nem.conf_initial_disclination_line()

nem.evolve_nematic_no_flow(5,method="ETD4RK")



for i in range(500):

    nem.evolve_nematic(10,method="ETD4RK")
    S,n = nem.calc_order_and_director()
    D = nem.calc_disclination_density_nematic()
    disclination_nodes = nem.calc_disclination_nodes_nem(charge_tolerance=0.0005)
    ax =nem.plot_field(S, vlim_symmetric=False,colormap='winter')
  #  nem.plot_nodes(disclination_nodes,ax=ax)
    plt.draw()
    plt.pause(0.01)
    cf.tool_save_plot_matplotlib(i)
    plt.clf()

cf.tool_make_animation_gif(i)

plt.show()