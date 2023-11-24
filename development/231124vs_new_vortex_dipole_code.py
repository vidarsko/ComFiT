import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=101,yRes=101)

theta = bec.calc_angle_field_vortex_dipole([bec.xmax/3,bec.ymax/3],[3*bec.xmax/3,3*bec.ymax/3])


bec.plot_angle_field(theta)
plt.show()