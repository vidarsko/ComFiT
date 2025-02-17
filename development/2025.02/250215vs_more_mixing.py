# import plotly.graph_objects as go

# fig = go.Figure()

# # First 2D plot in left half
# fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x', yaxis='y'))

# # Second 3D plot in right half
# fig.add_trace(go.Scatter3d(
#     x=[1, 2, 3], y=[1, 2, 3], z=[10, 20, 30],
#     mode='markers',
#     marker=dict(size=5)
# ))

# fig.update_layout(
#     xaxis=dict(domain=[0, 0.45]),
#     yaxis=dict(domain=[0, 1]),
#     scene=dict(domain=dict(x=[0.55, 1], y=[0, 1]))  # 3D plot domain
# )

# fig.show()


import plotly.graph_objects as go

fig = go.Figure()

# First plot in left half
fig.add_trace(go.Scatter(x=[1,2,3], y=[4,5,6], xaxis='x', yaxis='y'))

# Second plot in right half
fig.add_trace(go.Scatter(x=[1,2,3], y=[1,2,3], xaxis='x2', yaxis='y2'))

fig.update_layout(
    xaxis2=dict(domain=[0.55, 1]),
    yaxis2=dict(domain=[0, 1])
)

fig.show()