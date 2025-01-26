import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

sys2 = cf.BaseSystem(2, xlim=[0,2*np.pi], xRes=100, ylim=[0,2*np.pi], yRes=100)
f = np.sin(sys2.x + sys2.y)
f_f = sp.fft.fftn(f)
gradf = [sp.fft.ifftn(sys2.dif[0]*f_f), sp.fft.ifftn(sys2.dif[1]*f_f)]

sys2.plot_vector_field(gradf)

plt.show()