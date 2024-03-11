import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


bs = cf.BaseSystem(3, xRes=31, yRes=31, zRes=31)
field = bs.x

bs.plot_field(field, number_of_layers=3,alpha=1)
plt.show()