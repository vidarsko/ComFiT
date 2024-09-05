from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Create a 2x2 subplot, specifying 'domain' type for the pie chart subplot
fig = make_subplots(
    rows=2, 
    cols=2, 
    subplot_titles=("Plot 1", "Plot 2", "Plot 3", "Plot 4"),
    specs=[[{'type': 'xy'}, {'type': 'xy'}],  # First row: scatter/bar plots
           [{'type': 'domain'}, {'type': 'xy'}]]  # Second row: pie chart (domain type) and scatter plot
)

# Add traces to the subplots
fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], mode='lines', name='Line Plot'), row=1, col=1)
fig.add_trace(go.Bar(x=['A', 'B', 'C'], y=[7, 8, 9], name='Bar Plot'), row=1, col=2)
fig.add_trace(go.Pie(labels=['X', 'Y', 'Z'], values=[10, 15, 20], name='Pie Chart'), row=2, col=1)
fig.add_trace(go.Scatter(x=[1, 2, 3], y=[3, 1, 6], mode='markers', name='Scatter Plot'), row=2, col=2)

# Update layout
fig.update_layout(height=600, width=800, title_text="Subplots Example")

# Show the figure
fig.show()
