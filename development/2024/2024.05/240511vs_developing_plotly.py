import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Create a subplot grid of 2 rows by 1 column
fig = make_subplots(rows=2, cols=1)

# Add a scatter plot to the first subplot
fig.add_trace(
    go.Scatter(x=[1, 2, 3], y=[4, 5, 6]),
    row=1, col=1
)

# Add a bar plot to the second subplot
fig.add_trace(
    go.Bar(x=[1, 2, 3], y=[2, 3, 4]),
    row=2, col=1
)

# Update layout if needed (optional)
fig.update_layout(height=600, width=600, title_text="Subplots Example")
fig.show()

