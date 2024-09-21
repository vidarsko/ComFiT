import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

#2D system
bs = cf.BaseSystem(2,xRes=31,yRes=31)

# 1D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5)])
fig4 = bs.plot_vector_field(vector_field,spacing=3)

# 2D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5), bs.y*np.sin(bs.x/5)])
fig5 = bs.plot_vector_field(vector_field,spacing=5)


fig = cf.tool_make_subplots(1,2,fig4,fig5)

fig.show()