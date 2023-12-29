import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1)



nem.conf_initial_condition_disordered(noise_strength=.1)
#nem.evolve_nematic_no_flow(10)
nem.conf_active_channel(20,2)

plt.figure()
nem.plot_field(nem.alpha,cmap_symmetric=False,colormap="brg")
#plt.plot([nem.xmid,nem.xmid],[nem.ymin,nem.ymax-1],ls=':',color ='w')
#plt.plot([nem.xmid-10,nem.xmid-10],[nem.ymin,nem.ymax-1],ls=':',color ='w')
#plt.plot([nem.xmid+10,nem.xmid+10],[nem.ymin,nem.ymax-1],ls=':',color ='w')
plt.savefig("spatial_alpha.png",dpi=300)
plt.show()

#nem.evolve_nematic(400,"ETD4RK")
#D = nem.calc_defect_density_nematic()
#director =nem.calc_director()
#ax= nem.plot_field_velocity_and_director(D,nem.u,director)

#plt.show()