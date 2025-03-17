
import comfit as cf

bs = cf.BaseSystem(1, xlim=[-10,10])
field = bs.calc_Gaussian()

# fig, ax = bs.plot_field(field)
# bs.show(fig)

field_f = bs.fft(field)
fig, ax = bs.plot_complex_field(field_f, fourier=True)
bs.show(fig)