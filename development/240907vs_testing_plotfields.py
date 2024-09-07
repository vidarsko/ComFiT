import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from plotly.subplots import make_subplots


# fig = make_subplots(rows=1, cols=3, 
#                     subplot_titles=("1D", "2D", "3D"))

# 1D system

bs = cf.BaseSystem(1,xRes=31)
bs.plot_lib='plotly'
field = bs.x**2
fig=bs.plot_field(field)
fig.show()


# 2D system
bs = cf.BaseSystem(2,xRes=31,yRes=31)
bs.plot_lib='plotly'
field = bs.x**2+bs.y**2
fig = bs.plot_field(field)
fig.show()

# 3D system
bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
bs.plot_lib='plotly'
field = bs.x**2+bs.y**2+bs.z**2
fig = bs.plot_field(field, number_of_layers=3)
fig.show()