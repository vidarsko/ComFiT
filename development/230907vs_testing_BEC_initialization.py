import comfit as cf
import numpy as np

bec = cf.BEC(2,xRes=101,yRes=101)

bec.set_initial_condition_disordered()

angle=np.angle(bec.psi)

print(np.shape(angle))

print(bec.k[0] + bec.k[1] + bec.k[2])
bec.plot_angle_field(angle)
