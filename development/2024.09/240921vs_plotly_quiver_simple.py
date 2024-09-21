import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

bs = cf.BaseSystem(1, xRes=31)

# 2D vector field
vector_field = np.array([bs.x*np.cos(bs.x/5), bs.x*np.sin(bs.x/5)])
fig2 = bs.plot_vector_field(vector_field, spacing=5)

fig2.show()