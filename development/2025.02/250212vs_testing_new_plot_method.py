import comfit as cf
import numpy as np

cfi = cf.BoseEinsteinCondensate(2)
cfi.psi = cfi.x + 1j * cfi.y
cfi.plot_lib = 'plotly'

fig, axs = cfi.plot_subplots(2, 2, figsize=(10, 10))

cfi.plot_complex_field(cfi.psi, fig=fig, ax=axs[0][0])
# cfi.plot_field(abs(cfi.psi), fig=fig, ax=axs[0][1])
# cfi.plot_angle_field(np.angle(cfi.psi), fig=fig, ax=axs[1][0])
# cfi.show(fig)