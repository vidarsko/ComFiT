import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

sys3 = cf.BaseSystem(3,
                     xRes=21,dx=0.1, xmin=-1,
                     yRes=21,dy=0.1, ymin=-1,
                     zRes=21,dz=0.1, zmin=-1)
f = np.exp(-sys3.x**2 - sys3.y**2 - sys3.z**2)
sys3.plot_field(f,number_of_layers=3)
plt.show()