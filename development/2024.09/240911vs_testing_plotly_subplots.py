import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from plotly.subplots import make_subplots

fig = make_subplots(rows=2, cols=2, subplot_titles=('Plot 1', 'Plot 2', 'Plot 3', 'Plot 4'))

# fig = go.Figure()
# bs = cf.PhaseFieldCrystal1DPeriodic(10)
# field = np.sin(bs.x/bs.a0)

# fig = bs.plot_field(field, axis_equal=True)
fig.show() 

# fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]), 1, 1)
# fig.add_trace(go.Bar(x=[1, 2, 3], y=[6, 5, 4]), row=1, col=2)