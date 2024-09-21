import plotly.subplots
import plotly.graph_objs as go
import numpy as np

# Create the figure with subplots
fig = plotly.subplots.make_subplots(rows=2, cols=2)

# Upper Left plot (subplot at row=1, col=1)
fig.add_trace(
    go.Scatter(x=[1, 2, 3], y=[4, 5, 6], mode='lines', name='Plot 1'),
    row=1, col=1
)

# Upper Right plot (subplot at row=1, col=2)
fig.add_trace(
    go.Bar(x=[1, 2, 3], y=[6, 5, 4], name='Plot 2'),
    row=1, col=2
)

# Lower Left plot (subplot at row=2, col=1)
# Create data for the heatmap
z = np.random.rand(10, 10)

fig.add_trace(
    go.Heatmap(z=z, name='Heatmap'),
    row=2, col=1
)

# Set the aspect ratio to be equal for the lower left plot
fig.update_xaxes(scaleanchor="y3", scaleratio=1, row=2, col=1)
fig.update_yaxes(scaleanchor="x3", scaleratio=1, row=2, col=1)

# Lower Right plot (subplot at row=2, col=2)
fig.add_trace(
    go.Scatter(x=[1, 2, 3], y=[9, 8, 7], mode='lines', name='Plot 4'),
    row=2, col=2
)

# Update layout for better spacing
fig.update_layout(height=600, width=800, title_text="2x2 Subplots with Heatmap and Equal Axis")

# Display the figure
fig.show()
