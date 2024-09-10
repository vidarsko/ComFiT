import comfit as cf
import numpy as np

bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101)

bec.conf_initial_condition_disordered()

angle=np.angle(bec.psi)

print(np.shape(angle))

print(bec.calc_k2())
bec.plot_angle_field(angle)
