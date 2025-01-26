import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf

import numpy as np
from mayavi import mlab

bs = cf.BaseSystem(2,xRes=201,dx=0.05,yRes=201,dy=0.05)
z = np.sin(bs.x + 0*bs.y)

fig = mlab.figure(bgcolor=(1, 1, 1),fgcolor=(0.,0.,0.))
s = mlab.surf(bs.x, bs.y, z, colormap='viridis', figure=fig)

# Adding axes with labels
axes = mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', figure=fig, 
                 nb_labels=5, ranges=(0, 5, 0, 5, -1, 1))
axes.label_text_property.font_family = 'times'
axes.label_text_property.font_size = 10
# mlab.view(azimuth=-135, elevation=60,distance=15)
mlab.view(-135,60)

# Adding a colorbar
cb = mlab.colorbar(s, orientation='vertical', nb_labels=5, label_fmt='%.2f')
mlab.show()
