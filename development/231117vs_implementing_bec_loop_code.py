import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(3,xRes=21,yRes=21,zRes=21)
theta = bec.calc_angle_field_vortex_ring(normal_vector=[1,1,0])
#print(np.min(theta))
#print(np.max(theta))
#bec.plot_field(theta)
bec.plot_angle_field(theta)
plt.show()