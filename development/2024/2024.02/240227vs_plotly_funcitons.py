import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import plotly.graph_objects as go
import numpy as np

bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
bec.conf_initial_condition_disordered()
bec.evolve_relax(100)

x, y, z = np.meshgrid(bec.x, bec.y, bec.z)
values = np.array([0.5])

# Define the isosurface plot
fig = go.Figure(data=go.Isosurface(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value=values.flatten(),
    isomin=0,  # Minimum iso-value
    isomax=1,  # Maximum iso-value
    caps=dict(x_show=False, y_show=False)
))

fig.show()