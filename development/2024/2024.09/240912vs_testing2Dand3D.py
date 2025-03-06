
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

# fig = make_subplots(1,2, specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]])
# fig = make_subplots(1,2)
# fig = go.Figure()

# fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]),row=1, col=1)
# fig.add_trace(go.Scatter3d(x=[1, 2, 3], y=[4, 5, 6], z=[7, 8, 9]), row=1, col=1)
# fig.add_trace(go.Scatter3d(x=[3, 2, 1], y=[4, 5, 6], z=[7, 8, 9]))

# 3D scatter data
x_scatter = [1, 2, 3, 4, 5]
y_scatter = [5, 6, 2, 3, 1]
z_scatter = [2, 3, 3, 7, 8]

# Surface plot data
x = np.linspace(-5, 5, 50)
y = np.linspace(-5, 5, 50)
x, y = np.meshgrid(x, y)
z = np.sin(np.sqrt(x**2 + y**2))

# Add 3D scatter trace
# fig.add_trace(go.Scatter3d(x=x_scatter, y=y_scatter, z=z_scatter, 
#                            mode='markers', marker=dict(size=5, color='blue')))

# Add surface trace
# fig.add_trace(go.Surface(z=z, x=x, y=y, colorscale='viridis', opacity=0.7), row=1, col=1)
fig.add_trace(go.Surface(z=z, x=x, y=y, colorscale='viridis', opacity=0.7))


for trace in fig.data:
    # print('3d' in trace.type)
    print(trace.type)
fig.show()