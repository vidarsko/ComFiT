import comfit as cf
import random

plot_lib = 'matplotlib'

logmin = -15
logmax = 15

def make_plot(dim):
    if dim == 1:
        cfi = cf.BaseSystem(1,xRes=41)
        field = 10**random.randint(logmin,logmax)*cfi.x
    elif dim == 2:
        cfi = cf.BaseSystem(2,xRes=41,yRes=61)
        field = 10**random.randint(logmin,logmax)*(cfi.x + cfi.y)
    elif dim == 3:
        cfi = cf.BaseSystem(3,xRes=41,yRes=61,zRes=51)
        field = 10**random.randint(logmin,logmax)*(cfi.x + cfi.y + cfi.z)

    cfi.plot_lib = plot_lib

    return cfi, field

nrows = 2
ncols = 2
cfi = cf.BaseSystem(1,xRes=41, plot_lib=plot_lib)
fig, axs = cfi.plot_subplots(nrows, ncols)

for i in range(nrows):
    for j in range(ncols):
        cfi, field = make_plot(random.randint(3,3))
        cfi.plot_field_in_plane(field, fig=fig, ax=axs[i][j])

cfi.show(fig)

