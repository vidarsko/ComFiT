import comfit as cf
import numpy as np

bec = cf.BEC(2,xRes=101,yRes=101)

bec.set_initial_condition_disordered()

angle=np.angle(bec.psi)

print(np.shape(angle))
bec.plot_angle_field(angle)