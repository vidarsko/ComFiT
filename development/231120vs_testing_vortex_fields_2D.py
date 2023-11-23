import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=201,yRes=201)
#theta = bec.calc_angle_field_double_vortex(position1=[bec.xmax/3,bec.ymax/3],position2=[2*bec.xmax/3,bec.ymax/3])
#theta = theta + bec.calc_angle_field_double_vortex(position1=[2*bec.xmax/3,2*bec.ymax/3],position2=[bec.xmax/3,2*bec.ymax/3])
# theta = bec.calc_angle_field_double_vortex(position1=[bec.xmid-10,bec.ymid],position2=[bec.xmid+10,bec.ymid])
theta = bec.calc_angle_field_double_vortex()
bec.plot_angle_field(theta)
plt.show()