import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

N=6
res = 101
bec = cf.BoseEinsteinCondensate(2,xRes=(N+1)*res,yRes=(N+1)*res)

theta=0
# ny=3
# for nx in range(0,N+1):
#     print(nx,ny)
#     theta = theta + bec.calc_angle_field_vortex_dipole(position1=[res/3+nx*res,res/2+ny*res],position2=[2*res/3+nx*res,res/2+ny*res])


for ny in range(0, N + 1):
    for nx in range(0, N + 1):
        print(nx, ny)
        theta = theta + bec.calc_angle_field_vortex_dipole(position1=[res / 3 + nx * res, res / 2 + ny * res],
                                                   position2=[2 * res / 3 + nx * res, res / 2 + ny * res])

#print(theta)
bec.plot_angle_field(theta)
plt.show()