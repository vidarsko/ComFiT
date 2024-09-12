from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
fig = make_subplots(1, 2, horizontal_spacing=0.15)
fig.add_trace(go.Heatmap(z=np.random.random((5, 5)), colorbar_x=0.45), 1, 1)
fig.add_trace(go.Heatmap(z=np.random.random((5, 5)) + 3), 1, 2)
fig.show()