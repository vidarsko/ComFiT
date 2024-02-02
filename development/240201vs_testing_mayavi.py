import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
from mayavi import mlab

bs = cf.BaseSystem(2, xRes=201, dx=0.05, yRes=201, dy=0.05)

z = np.sin(bs.x + bs.y)

fig = mlab.figure(bgcolor=(1, 1, 1),fgcolor=(0.,0.,0.))
s = mlab.surf(bs.x, bs.y, z, colormap='viridis', figure=fig)

# Adding axes with labels
axes = mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', figure=fig)
axes.label_text_property.font_family = 'times'
axes.label_text_property.font_size = 10



mlab.show()
