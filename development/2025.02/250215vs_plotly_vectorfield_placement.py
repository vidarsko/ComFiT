import comfit as cf
import random
import numpy as np

plot_lib = 'plotly'

logmin = -15
logmax = 15

def make_plot(dim,vector_dim):
    if dim == 1:
        cfi = cf.BaseSystem(1,xRes=41)
        field = np.zeros((vector_dim,cfi.xRes))
        order_of_magnitude = random.randint(logmin,logmax)
        for i in range(vector_dim):
            field[i] = 10**order_of_magnitude*np.sin(cfi.x/cfi.xmax)
    elif dim == 2:
        cfi = cf.BaseSystem(2,xRes=41,yRes=61)
        field = np.zeros((vector_dim,cfi.xRes,cfi.yRes))
        order_of_magnitude = random.randint(logmin,logmax)
        for i in range(vector_dim):
            field[i] = 10**order_of_magnitude*(cfi.x/cfi.xmax + np.sin(cfi.y/cfi.ymax))
    elif dim == 3:
        cfi = cf.BaseSystem(3,xRes=41,yRes=61,zRes=51)
        field = np.zeros((vector_dim,cfi.xRes,cfi.yRes,cfi.zRes))
        order_of_magnitude = random.randint(logmin,logmax)
        for i in range(vector_dim):
            field[i] = 10**order_of_magnitude*(cfi.x/cfi.xmax + np.cos(cfi.y/cfi.ymax) + np.sin(cfi.z/(2*cfi.zmax)))

    cfi.plot_lib = plot_lib

    return cfi, field

nrows = 3
ncols = 3
cfi = cf.BaseSystem(1,xRes=41, plot_lib=plot_lib)
fig, axs = cfi.plot_subplots(nrows, ncols)
# print(axs)


for i in range(nrows):
    for j in range(ncols):
        cfi, field = make_plot(i+1, j+1)
        cfi.plot_vector_field(field, fig=fig, ax=axs[i][j])

# print(axs)
cfi.show(fig)

