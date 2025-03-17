import numpy as np
import scipy.fft as sp_fft

# Example field (replace with your actual data)
N = 128  # number of grid points
dx = 0.1  # spatial step size
x0 = 5.0  # physical offset of the field (starting point)
x = np.arange(N) * dx + x0  # physical coordinates
field = np.sin(2 * np.pi * x / 10)  # example field with offset

# Compute the FFT
fourier_field = sp_fft.fftn(field)

# Generate wave numbers (assuming periodic domain of length L = N * dx)
k = sp_fft.fftfreq(N, d=dx)  # wave numbers in cycles per unit length
k = 2 * np.pi * k  # convert to angular frequency (radians per unit length)

# Apply phase correction for the offset x0
phase_correction = np.exp(1j * k * x0)  # e^(i k x0)
corrected_fourier_field = fourier_field * phase_correction

# Now you can proceed with plotting or further analysis

