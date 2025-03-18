
import comfit as cf
import numpy as np

bs = cf.BaseSystem(1, xlim=[-10,10], xRes=101)

field = 0.5*bs.calc_Gaussian(position=1)
field_f = bs.fft(field)

fig, ax = bs.plot_subplots(1,2)
bs.plot_field(field, ax=ax[0], fig=fig, title='Real field')
bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig, title='Fourier transform')

bs.plot_save(fig)

# bs.plot_save(fig)
# bs = cf.BaseSystem(2, xlim=[-10,10], ylim=[-10,10])
# field = bs.calc_Gaussian(position=(0,0))

# fig, ax = bs.plot_subplots(2,1)

# bs.plot_field(field, ax=ax[0], fig=fig)

# field_f = bs.fft(field)
# fig, ax = bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig)

# bs.show(fig)


# bs = cf.BaseSystem(3, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
# field = bs.calc_Gaussian(position=(0,0,0))

# fig, ax = bs.plot_subplots(2,1)

# bs.plot_field(field, ax=ax[0], fig=fig)

# field_f = bs.fft(field)
# fig, ax = bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig)

# bs.show(fig)
