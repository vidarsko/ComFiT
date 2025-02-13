import comfit as cf
import numpy as np

cfi = cf.BoseEinsteinCondensate(2)
cfi.psi = cfi.x + 1j * cfi.y

cfi.plot_lib = 'plotly'

nrows = 2
ncols = 3

fig, axs = cfi.plot_subplots(nrows, ncols, figsize=(10, 10))

for i in range(nrows):
    for j in range(ncols):
        cfi.plot_complex_field(cfi.psi, fig=fig, ax=axs[i][j])
# cfi.plot_complex_field(cfi.psi, fig=fig, ax=axs[0][0])
# fig, ax = cfi.plot_complex_field(cfi.psi)
# cfi.plot_field(abs(cfi.psi), fig=fig, ax=axs[0][1])
# cfi.plot_angle_field(np.angle(cfi.psi), fig=fig, ax=axs[1][0])
cfi.show(fig)