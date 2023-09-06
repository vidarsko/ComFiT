from comfit.core.base_system import BaseSystem

bec = BaseSystem(2, xRes=101, yRes=101)

field = bec.calc_angle_field_double_vortex()

bec.plot_angle_field(field)
