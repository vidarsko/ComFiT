import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101)

#Figure a, just adding
#theta = bec.calc_angle_field_vortex_dipole([0.45*bec.xmax,0],[0.75*bec.xmax,bec.ymid])
# theta = 0
# theta = theta + bec.calc_angle_field_single_vortex([0.5*bec.xmax,0.5*bec.ymid], charge=-1)
# theta = theta + bec.calc_angle_field_single_vortex([0.98*bec.xmax,0.5*bec.ymid], charge=1)
# theta = np.mod(theta+np.pi, 2 * np.pi) - np.pi

theta = bec.calc_angle_field_vortex_dipole([0.48*bec.xmax,0],[0.74*bec.xmax,0.5*bec.ymid])

ax = bec.plot_angle_field(theta,colorbar=False)
plt.gcf().set_size_inches(2.5,2.5)
plt.gcf().tight_layout()
plt.savefig(r'C:\Users\vidar\UiO Dropbox\Vidar Skogvoll\Apps\Overleaf\230531 - Code documentation\Figures\NumericalImplementationOfAngleFields\test.png',dpi=300)