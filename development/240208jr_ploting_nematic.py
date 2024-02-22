import numpy as np
import matplotlib.pyplot as plt
import comfit as cf
from mayavi import mlab

nem = cf.NematicLiquidCrystal(3,xRes=31,yRes=31,zRes = 31,dx=1,dy=1,dt=0.1,alpha=-1.5,C=1)



nem.conf_insert_disclination_line(angle = 1)

nem.evolve_nematic_no_flow(50,method="ETD4RK")


#nem.evolve_nematic(800,method="ETD4RK")

S,n = nem.calc_order_and_director()
nem.plot_nematic_3D(S,director=True)
mlab.show()