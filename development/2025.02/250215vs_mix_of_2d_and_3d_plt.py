import comfit as cf
import numpy as np

bec = cf.BaseSystem(1, xRes=20)

bec.plot_lib = 'plotly'

field = (np.random.rand(bec.xRes))

fig, axs = bec.plot_subplots(1,2)

bec.plot_field(field, fig=fig, ax=axs[1])

bec2 = cf.BaseSystem(3, xRes=10, yRes=10, zRes=10)

bec2.plot_lib = 'plotly'

# field = (np.random.rand(bec2.xRes, bec2.yRes, bec2.zRes))

# bec2.plot_field(field, fig=fig, ax=axs[0])
bec.plot_field(field, fig=fig, ax=axs[0])

bec.show(fig)