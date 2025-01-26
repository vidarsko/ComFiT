import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

bs = cf.BaseSystem(2)
field = np.zeros((bs.xRes, bs.yRes))

field[:,:] = np.nan
# field[0,0] = 1
# field[50,50] = 0

fig = bs.plot_field(field)
cf.tool_save_plot(1,fig)
# fig.show()

