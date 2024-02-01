from mayavi import mlab
import numpy as np

x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
x, y = np.meshgrid(x, y)
z = np.sin(x**2 + y**2)

# Create a 3D surface plot
mlab.figure(bgcolor=(1, 1, 1))
surf = mlab.surf(x, y, z)

# Close the figure
mlab.show()
