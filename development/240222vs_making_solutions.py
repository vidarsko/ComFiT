import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bs = cf.BaseSystem(1, xRes=12)

print(bs.x)
print(bs.xmidi)
print(bs.xmax)

# f = np.sin(bs.x+bs.y)
# bs.plot_field(f)
# plt.show()