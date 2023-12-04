import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

plt.figure(figsize=(5.52, 5))

mappable = plt.cm.ScalarMappable(cmap=cf.tool_colormap_angle())
mappable.set_array([])
mappable.set_clim(-np.pi, np.pi)
cbar = plt.colorbar(mappable,orientation='horizontal')
cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

plt.savefig('231204vs_color_scheme.png',dpi=300)