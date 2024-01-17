import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(3,xRes=31,yRes=31,zRes=31,dx=1,dy=1,dt=0.1,alpha=-1.5)



nem.conf_initial_condition_disordered(noise_strength=1.0)
nem.evolve_nematic_no_flow(10,method="ETD4RK")

nem.evolve_nematic(2000,method="ETD4RK")

S_0 =nem.calc_S()
S,n = nem.calc_order_and_director()




plt.figure()
nem.plot_field(S ,cmap_symmetric=False,colormap="brg",number_of_layers=2)
#nem.plot_nematic_3D(n,S)
#plt.plot([nem.xmid,nem.xmid],[nem.ymin,nem.ymax-1],ls=':',color ='w')
#plt.plot([nem.xmid-10,nem.xmid-10],[nem.ymin,nem.ymax-1],ls=':',color ='w')
#plt.plot([nem.xmid+10,nem.xmid+10],[nem.ymin,nem.ymax-1],ls=':',color ='w')
#plt.savefig("3D_nematic_test.png",dpi=300)
plt.show()

#nem.evolve_nematic(400,"ETD4RK")
#D = nem.calc_defect_density_nematic()
#director =nem.calc_director()
#ax= nem.plot_field_velocity_and_director(D,nem.u,director)

#plt.show()