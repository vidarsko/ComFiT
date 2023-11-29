import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=31,yRes=31,zRes=31)
#bec.set_initial_condition_vortex_ring(normal_vector=[1,1,0])
#bec.evolve_relax_BEC(200)
#bec.plot_field(abs(bec.psi))
theta = bec.calc_angle_field_vortex_ring(normal_vector=[1,0,1])
bec.plot_angle_field(theta)
plt.show()