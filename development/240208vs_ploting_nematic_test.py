import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
from mayavi import mlab

nem = cf.NematicLiquidCrystal(3,xRes=31,yRes=31,zRes = 31,dx=1,dy=1,dt=0.1,alpha=-1.5,C=1)



nem.conf_initial_condition_disordered(noise_strength=1.0)
nem.evolve_nematic_no_flow(5,method="ETD4RK")


nem.evolve_nematic(10,method="ETD4RK")

S,n = nem.calc_order_and_director()
nem.plot_nematic_3D(S,director=True, Flow=True)
mlab.show()