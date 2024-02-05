import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt


#1D
# Plot normal field
# bs = cf.BaseSystem(1, xRes = 101, dx=0.1)
# y = np.sin(bs.x)

# bs.plot_field(y)

#2D
# Plot normal field
bs = cf.BaseSystem(2,xRes = 101, dx = 0.1, yRes = 101, dy = 0.1)
field = np.sin(bs.x) * np.sin(bs.y)

bs.plot_field(field, xlim=[1,5], ylim=[1,2])
plt.show()

#3D
# Plot normal field
# bs = cf.BaseSystem(3,xRes = 101, dx = 0.1, yRes = 101, dy = 0.1, zRes = 101, dz = 0.1)
# field = np.sin(bs.x) * np.sin(bs.y) * np.sin(bs.z)

# bs.plot_field(field,'suptite')
# plt.show()