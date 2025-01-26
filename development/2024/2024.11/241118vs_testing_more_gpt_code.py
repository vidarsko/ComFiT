
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

class SolveStaticField2D(cf.BaseSystem):
    def __init__(self, dim, **kwargs):
        super().__init__(dim, **kwargs)
        
    def solve_field(self, V):
        # Broadcast x to 2D (x.shape = (xRes, 1), y.shape = (1, yRes))
        psi = np.zeros_like(V)  # Initialize field
        non_zero_indices = self.x != 0  # Only divide where x != 0
        psi[non_zero_indices, :] = V[non_zero_indices, :] / (2 * self.x[non_zero_indices, :])
        return psi

# Define spatial domain and resolution
dim = 2
xmin, xmax = -10, 10
ymin, ymax = -10, 10
xRes, yRes = 100, 100

# Create instance of the 2D system
system = SolveStaticField2D(dim=dim, xmin=xmin, xmax=xmax, xRes=xRes,
                            ymin=ymin, ymax=ymax, yRes=yRes)

# Define V(x, y) as a 2D Gaussian potential
V = np.exp(-(system.x**2 + system.y**2))  # 2D Gaussian centered at (0, 0)

# Solve for psi(x, y)
psi = system.solve_field(V)

# Plot the result using matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Plot V(x, y)
c1 = ax[0].contourf(system.x + system.y, V, levels=50, cmap='viridis')
ax[0].set_title("V(x, y)")
fig.colorbar(c1, ax=ax[0])

# Plot psi(x, y)
c2 = ax[1].contourf(system.x + system.y, psi, levels=50, cmap='plasma')
ax[1].set_title("psi(x, y)")
fig.colorbar(c2, ax=ax[1])

plt.show()
