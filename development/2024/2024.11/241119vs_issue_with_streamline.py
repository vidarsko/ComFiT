import numpy as np
import plotly.figure_factory as ff

# Define the grid
x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
X,Y = np.meshgrid(x,y)

# Define the vector field
U = np.cos(X) * Y
V = np.sin(Y) * X

# Create streamline plot
fig = ff.create_streamline(x,y, U, V)

fig.write_image("streamline_plot.png")
