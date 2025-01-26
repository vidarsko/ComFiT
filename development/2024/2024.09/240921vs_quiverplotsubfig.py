import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

#1D system
bs = cf.BaseSystem(1,xRes=31)

# 1D vector field
vector_field = np.array([bs.x*np.cos(bs.x/5)])
fig1 = bs.plot_vector_field(vector_field, spacing=1)

# 2D vector field
vector_field = np.array([bs.x*np.cos(bs.x/5), bs.x*np.sin(bs.x/5)])
fig2 = bs.plot_vector_field(vector_field, spacing=2)

# 3D vector field
vector_field = np.array([bs.x*np.cos(bs.x/5), bs.x*np.sin(bs.x/5), bs.x*np.cos(bs.x/5)])
fig3 = bs.plot_vector_field(vector_field,spacing=3)

#2D system
bs = cf.BaseSystem(2,xRes=31,yRes=31)

# 1D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5)])
fig4 = bs.plot_vector_field(vector_field,spacing=3)

# 2D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5), bs.y*np.sin(bs.x/5)])
fig5 = bs.plot_vector_field(vector_field,spacing=5)

# 3D vector field
vector_field = np.array([bs.x*np.cos(bs.y/5), bs.y*np.sin(bs.x/5), bs.x*np.cos(bs.y/5)])
fig6 = bs.plot_vector_field(vector_field, spacing=3)

# 3D system
bs = cf.BaseSystem(3,xRes=11,yRes=11,zRes=11)

# 1D vector field
vector_field = np.array([bs.z+bs.x*np.cos(bs.y/5)])
fig7 = bs.plot_vector_field(vector_field,spacing=3)

# 2D vector field
vector_field = np.array([bs.z+ bs.x*np.cos(bs.y/5), bs.z + bs.y*np.sin(bs.x/5)])
fig8 = bs.plot_vector_field(vector_field,spacing=5)

# 3D vector field
vector_field = np.array([bs.z+ bs.x*np.cos(bs.y/5), bs.z + bs.y*np.sin(bs.x/5), -bs.z + bs.x*np.cos(bs.y/5)])
fig9 = bs.plot_vector_field(vector_field,spacing=3)


fig = cf.tool_make_subplots(3,3,fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9)
fig.show()


