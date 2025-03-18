
import comfit as cf

# bs = cf.BaseSystem(1, xlim=[-10,10], xRes=1001)
# # bs.plot_lib = 'matplotlib'
# field = bs.calc_Gaussian(position=0)

# fig, ax = bs.plot_subplots(2,1)

# field_f = bs.fft(field)

# bs.plot_field(field, ax=ax[0], fig=fig)


# bs.plot_field(field, ax=ax[1], fig=fig)

# bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig)




# bs.plot_save(fig)
bs = cf.BaseSystem(2, xlim=[-10,10], ylim=[-10,10])
field = bs.calc_Gaussian(position=(0,0))

fig, ax = bs.plot_subplots(2,1)

bs.plot_field(field, ax=ax[0], fig=fig)

field_f = bs.fft(field)
fig, ax = bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig)

bs.show(fig)


# bs = cf.BaseSystem(3, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
# field = bs.calc_Gaussian(position=(0,0,0))

# fig, ax = bs.plot_subplots(2,1)

# bs.plot_field(field, ax=ax[0], fig=fig)

# field_f = bs.fft(field)
# fig, ax = bs.plot_complex_field(field_f, fourier=True, ax=ax[1], fig=fig)

# bs.show(fig)
