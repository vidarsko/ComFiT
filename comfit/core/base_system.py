import numpy as np

class BaseSystem:
    def __init__(self, dimension, x_resolution, dx=1, dy=1, dz=1, dt=0.1):
        self.dim = dimension
        self.xRes = x_resolution
        self.yRes = 1
        self.zRes = 1
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.dt = dt

        if self.dim > 1:
            self.yRes = dy
        if self.dim > 2:
            self.zRes = dz

        if self.dim not in [1, 2, 3]:
            raise ValueError('Dimension must be 1, 2, or 3.')

        self.x = np.arange(0, self.xRes * self.dx, self.dx)
        self.y = np.array([0])
        self.z = np.array([0])

        if self.dim > 1:
            self.y = np.arange(0, self.yRes * self.dy, self.dy)
        if self.dim > 2:
            self.z = np.arange(0, self.zRes * self.dz, self.dz)

        self.Res = self.xRes * self.yRes * self.zRes
        self.dims = [self.yRes, self.xRes, self.zRes]

        self.a0 = 1  # System length PFCscale, set to 1 unless changed later

        # Helpful midpoints and their indices
        self.xmidi = (1 + self.xRes) // 2
        self.xmid = self.x[self.xmidi]

        self.ymidi = (1 + self.yRes) // 2
        self.ymid = self.y[self.ymidi]

        self.zmidi = (1 + self.zRes) // 2
        self.zmid = self.z[self.zmidi]

        self.midi = self.xRes * self.yRes * (self.zmidi - 1) + self.yRes * (self.xmidi - 1) + self.ymidi
        self.rmid = [self.xmid, self.ymid, self.zmid]

        # Max positions
        self.xmax = self.x[-1] + self.dx
        self.ymax = self.y[-1] + self.dy if self.dim > 1 else 0
        self.zmax = self.z[-1] + self.dz if self.dim > 2 else 0

        # Fourier modes
        self.ki = [self.calc_wavenums(self.x), None, None]
        if self.dim > 1:
            self.ki[1] = self.calc_wavenums(self.y)
        if self.dim > 2:
            self.ki[2] = self.calc_wavenums(self.z)

        # Derivatives
        self.dif = [1j * ki for ki in self.ki]

        self.dV = self.dx
        if self.dim > 1:
            self.dV *= self.dy
        if self.dim > 2:
            self.dV *= self.dz

        self.xmin = 0
        self.ymin = 0
        self.zmin = 0

        self.rmin = [self.xmin, self.ymin, self.zmin]
        self.rmax = [self.xmax, self.ymax, self.zmax]

    @staticmethod
    def calc_wavenums(x):
        """
        Calculates the wavenumbers corresponding to the input position vectors given by x.

        Parameters:
        - x : numpy array
            1D array of x-positions.

        Returns:
        - k : numpy array
            1D array of wavenumbers with all the modes for the given x-array,
            assuming periodicity from x[0] to x[0] over n intervals.

        Example:
        x = np.array([-10, -5, 0, 5, 10])
        k = instance_of_BaseSystem.calc_wavenums(x)
        print(k)
        # Output: [ 0.          0.25132741  0.50265482 -0.50265482 -0.25132741]
        """
        n = len(x)

        high = n // 2
        low = - (n - 1) // 2

        L = n * (x[1] - x[0])

        k = np.concatenate((np.arange(0, high+1), np.arange(low, 0))) * 2 * np.pi / L

        if n % 2 == 0:
            k[n//2] = -k[n//2]

        return k



# Usage:
# sys = BaseSystem(2, 100, dy=2)
