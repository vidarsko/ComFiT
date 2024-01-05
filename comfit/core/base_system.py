import numpy as np
import matplotlib.pyplot as plt
from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_parula
from comfit.tools.tool_create_orthonormal_triad import tool_create_orthonormal_triad
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting
from skimage.measure import marching_cubes
from matplotlib.tri import Triangulation
import matplotlib.colors as mcolors
import scipy as sp

class BaseSystem:
    def __init__(self, dimension,
                 xRes=101, dx=1.0, xmin=0,
                 yRes=101, dy=1.0, ymin=0,
                 zRes=101, dz=1.0, zmin=0,
                 dt=0.1, **kwargs):
        """
        Initialize the class with the given parameters.
        Parameters:
            dimension (int): The dimension of the system. Must be 1, 2, or 3.
            xRes (int): The resolution of the x-axis.
            dx (float): The spacing between points on the x-axis.
            yRes (int): The resolution of the y-axis.
            dy (float): The spacing between points on the y-axis.
            zRes (int): The resolution of the z-axis.
            dz (float): The spacing between points on the z-axis.
            dt (float): The time step.
            **kwargs: Additional keyword arguments.
        Raises:
            ValueError: If the dimension is not 1, 2, or 3.
        Returns:
            None
        """
        self.dim = dimension
        self.xRes = xRes
        self.yRes = 1
        self.zRes = 1
        self.time = 0
        if dimension > 1:
            self.yRes = yRes

        if dimension > 2:
            self.zRes = zRes

        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.dt = dt

        self.xmin = xmin
        self.xmax = xmin + self.xRes * self.dx

        self.ymin = ymin
        self.ymax = ymin + self.yRes * self.dy

        self.zmin = zmin
        self.zmax = zmin + self.zRes* self.dz

        if self.dim not in [1, 2, 3]:
            raise ValueError('Dimension must be 1, 2, or 3.')

        self.x = np.linspace(self.xmin, self.xmax-self.dx, self.xRes)
        self.y = np.linspace(self.ymin, self.ymax-self.dy, self.yRes)
        self.z = np.linspace(self.zmin, self.zmax-self.dz, self.zRes)

        self.Res = self.xRes * self.yRes * self.zRes
        if self.dim == 1:
            self.dims = 1
        elif self.dim == 2:
            self.dims = [self.xRes, self.yRes]
        elif self.dim == 3:
            self.dims = [self.xRes, self.yRes, self.zRes]

        self.a0 = 1  # System length scale

        # Helpful midpoints and their indices
        self.xmidi = (1 + self.xRes) // 2 - 1
        self.xmid = self.x[self.xmidi]

        self.ymidi = (1 + self.yRes) // 2 - 1
        self.ymid = self.y[self.ymidi]

        self.zmidi = (1 + self.zRes) // 2 - 1
        self.zmid = self.z[self.zmidi]

        self.midi = self.xRes * self.yRes * (self.zmidi - 1) + self.yRes * (self.xmidi - 1) + self.ymidi
        if self.dim == 1:
            self.rmid = self.xmid
        elif self.dim == 2:
            self.rmid = [self.xmid, self.ymid]
        elif self.dim == 3:
            self.rmid = [self.xmid, self.ymid, self.zmid]

        # Fourier modes
        self.k = [self.calc_wavenums(self.x)]
        if self.dim == 2:
            self.k[0] = self.k[0].reshape(self.xRes, 1)
            self.k.append(self.calc_wavenums(self.y).reshape(1, self.yRes))
        elif self.dim == 3:
            self.k[0] = self.k[0].reshape(self.xRes, 1, 1)
            self.k.append(self.calc_wavenums(self.y).reshape(1, self.yRes, 1))
            self.k.append(self.calc_wavenums(self.z).reshape(1, 1, self.zRes))

        # Reshape the position vectors
        if self.dim == 2:
            self.x = self.x.reshape((self.xRes, 1))
            self.y = self.y.reshape((1, self.yRes))
        elif self.dim == 3:
            self.x = self.x.reshape((self.xRes, 1, 1))
            self.y = self.y.reshape((1, self.yRes, 1))
            self.z = self.z.reshape((1, 1, self.zRes))


        # Derivatives
        self.dif = [1j * ki for ki in self.k]

        self.dV = self.dx
        if self.dim > 1:
            self.dV *= self.dy
        if self.dim > 2:
            self.dV *= self.dz

        self.rmin = [self.xmin, self.ymin, self.zmin]
        self.rmax = [self.xmax, self.ymax, self.zmax]

    # CALCULATION FUNCTIONS

    # Calculation of angle fields for vortices of different types
    def calc_angle_field_single_vortex(self,
                                       position=None,
                                       charge=1):
        """
        Calculate the angle field due to a single vortex.
        Args:
            position (list, optional): The position of the vortex. Defaults to None.
            charge (int, optional): The charge of the vortex. Defaults to 1.
        Raises:
            Exception: If the dimension of the system is not 2.
        Returns:
            numpy.ndarray: The angle field calculated for the vortex.
        """
        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]

        x = self.x.reshape((self.xRes, 1))
        y = self.y.reshape((1, self.yRes))

        theta = charge * np.arctan2(y - position[1], x - position[0])

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_angle_field_vortex_dipole(self,
                                       dipole_vector=None,
                                       dipole_position=None):
        """
        Calculates the angle field for a double vortex system.
        Args:
            dipole_vector (list, optional): The dipole vector. Defaults to None.
            dipole_position (list, optional): The position the center of mass of the dipole. Defaults to None.
        Raises:
            Exception: If the dimension of the system is not 2.
        Returns:
            float: The calculated angle field for the double vortex system.
        """

        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if dipole_vector is None:
            dipole_vector = [self.xmax / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        # Add the vortices to the theta-field
        print("dipole position", dipole_position)
        print("dipole vector", dipole_vector)
        
        theta = np.zeros(self.dims)
        theta += self.calc_angle_field_single_vortex(dipole_position - np.array(dipole_vector) / 2,
                                                            charge=-1)
        theta += self.calc_angle_field_single_vortex(dipole_position + np.array(dipole_vector) / 2, 
                                                            charge=1)

        # Convert the field to a complex field to make it fit the periodic boundary conditions
        amp = np.exp(1j * theta)

        # Filter the angle field
        width = 0.2 * np.min([self.xmax, self.ymax])
        radius = 0.4 * np.min([self.xmax, self.ymax])

        r2 = (self.x.reshape((self.xRes, 1)) - self.xmid) ** 2 + (self.y.reshape((1, self.yRes)) - self.ymid) ** 2
        filter = (1 + np.tanh((radius ** 2 - r2) / width ** 2)) / 2
        amp = amp * filter + (1 - filter)

        theta = np.angle(amp)

        # Roll the field so that the dipole center is at the desired position.
        Rx = round((self.rmid[0] - dipole_position[0]) / self.dx)
        theta = np.roll(theta, -Rx, axis=0)
        Ry = round((self.rmid[1] - dipole_position[1]) / self.dy)
        theta = np.roll(theta, -Ry, axis=1)

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_angle_field_vortex_ring(self, position=None, radius=None, normal_vector=[0, 0, 1]):

        if position is None:
            position = self.rmid

        if radius is None:
            radius = np.min([self.xmax, self.ymax, self.zmax]) / 3

        if radius > np.min([self.xmax, self.ymax, self.zmax]) / 3:
            print("Warning: The radius of the suggested vortex ring is large."
                  "This can cause unwanted boundary effects.")

        n = normal_vector / np.linalg.norm(np.array(normal_vector))
        [X, Y, Z] = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Add the vortices to the theta-field
        theta = 0
        position = np.array(position)

        position = self.rmid

        m2 = n[0] * (X - position[0]) \
             + n[1] * (Y - position[1]) \
             + n[2] * (Z - position[2])

        m1 = np.sqrt(
            (X - position[0] - m2 * n[0]) ** 2
            + (Y - position[1] - m2 * n[1]) ** 2
            + (Z - position[2] - m2 * n[2]) ** 2
        )

        theta = theta + np.arctan2(m2, m1 + radius)
        theta = theta + np.arctan2(m2, m1 - radius)

        # Convert the field to a complex field to make it fit the periodic boundary conditions
        amp = np.exp(1j * theta)

        # Filter the angle field
        width = 0.2 * np.min([self.xmax, self.ymax, self.zmax])
        radius = 0.4 * np.min([self.xmax, self.ymax, self.zmax])
        # TODO: This radius shares name with the one defining the ring. May cause trouble down the line (Vidar 01.12.23)

        r2 = (self.x.reshape((self.xRes, 1, 1)) - self.xmid) ** 2 \
             + (self.y.reshape((1, self.yRes, 1)) - self.ymid) ** 2 \
             + (self.z.reshape((1, 1, self.zRes)) - self.zmid) ** 2

        filter = (1 + np.tanh((radius ** 2 - r2) / width ** 2)) / 2  # Goes from 1 to 0 continuously
        amp = amp * filter + (1 - filter)

        theta = np.angle(amp)

        # Roll the field so that the ring center is at the desired position.
        Rx = round((self.rmid[0] - position[0]) / self.dx)
        theta = np.roll(theta, -Rx, axis=0)
        Ry = round((self.rmid[1] - position[1]) / self.dy)
        theta = np.roll(theta, -Ry, axis=1)
        Rz = round((self.rmid[2] - position[2]) / self.dz)
        theta = np.roll(theta, -Rz, axis=2)

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_wavenums(self, x):
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
        k = instance_of_BaseSystem.calc_wavenums(self,x)
        print(k)
        # Output: [ 0.          0.25132741  0.50265482 -0.50265482 -0.25132741]
        """
        n = len(x)

        high = (n - 1) // 2
        low = - (n // 2)

        l = n * (x[1] - x[0])

        k = np.concatenate((np.arange(0, high + 1), np.arange(low, 0))) * 2 * np.pi / l

        return k

    def calc_k2(self):
        return sum([self.k[i] ** 2 for i in range(len(self.k))])

    def calc_gaussfilter_f(self, a0=None):

        if a0 is None:
            a0 = self.a0

        return np.exp(-1 / 2 * a0 ** 2 * self.calc_k2())

    def calc_determinant_field(self, psi):
        """
        Calculate the determinant transformation of a given field

        Parameters:
            psi (list): A list of two psi fields.

        Returns:
            float: The defect density of the psi field.
        """
        
        if self.dim == 2:
            if len(psi) == 2:
                psi_f = [sp.fft.fftn(psi[0]), sp.fft.fftn(psi[1])]

                return np.real(
                    sp.fft.ifftn(self.dif[0] * psi_f[0]) * sp.fft.ifftn(self.dif[1] * psi_f[1]) -
                    sp.fft.ifftn(self.dif[1] * psi_f[0]) * sp.fft.ifftn(self.dif[0] * psi_f[1]))
            
        elif self.dim == 3:
            if len(psi) == 2:
                psi_f = [sp.fft.fftn(psi[0]), sp.fft.fftn(psi[1])]

                result = np.array([
                    sp.fft.ifftn(self.dif[1] * psi_f[0]) * sp.fft.ifftn(self.dif[2] * psi_f[1]) -
                    sp.fft.ifftn(self.dif[2] * psi_f[0]) * sp.fft.ifftn(self.dif[1] * psi_f[1]),
                    -sp.fft.ifftn(self.dif[0] * psi_f[0]) * sp.fft.ifftn(self.dif[2] * psi_f[1]) +
                    sp.fft.ifftn(self.dif[2] * psi_f[0]) * sp.fft.ifftn(self.dif[0] * psi_f[1]),
                    sp.fft.ifftn(self.dif[0] * psi_f[0]) * sp.fft.ifftn(self.dif[1] * psi_f[1]) -
                    sp.fft.ifftn(self.dif[1] * psi_f[0]) * sp.fft.ifftn(self.dif[0] * psi_f[1])
                ], dtype='float64')
                #TODO: verify that this is correct to specify as float64. (vidar 26.12.23)

                return np.array(result)

    def calc_defect_density(self, psi, psi0=1):
        """
        Calculate the defect density of a given psi field.

        Parameters:
            psi (list): A list of two psi fields.
            psi0 (float, optional): The value of psi_0. Defaults to 1.

        Returns:
            float: The defect density of the psi field.
        """

        return 1 / (np.pi * psi0 ** 2) * self.calc_determinant_field(psi)

    def calc_defect_density_singular(self, psi, psi0=1):
        """
        Calculate the singular defect density for a given psi field.

        Parameters:
            psi (float): The value of psi.
            psi0 (float, optional): The reference value of psi. Defaults to 1.
        Returns:
            float: The defect density for the given psi value.
        """
        return self.calc_defect_density(psi, 1) * self.calc_delta_function(psi, psi0)

    def calc_defect_velocity_field(self, psi, dt_psi):

        if self.dim == 2:
            if len(psi) == 2:
                # Parameters to exclude region
                threshold = 0.4

                psi_f = [sp.fft.fftn(psi[0]), sp.fft.fftn(psi[1])]

                dx_psi0 = sp.fft.ifftn(self.dif[0] * psi_f[0])
                dy_psi1 = sp.fft.ifftn(self.dif[1] * psi_f[1])
                dx_psi1 = sp.fft.ifftn(self.dif[0] * psi_f[1])
                dy_psi0 = sp.fft.ifftn(self.dif[1] * psi_f[0])

                denominator = np.real(dx_psi0 * dy_psi1 - dx_psi1 * dy_psi0)

                denominator_max = np.max(abs(denominator))
                region_to_ignore = abs(denominator) < threshold * denominator_max

                Vx = -2 * np.real(dt_psi[0] * dy_psi1 - dt_psi[1] * dy_psi0) / denominator
                Vy = -2 * np.real(-dt_psi[0] * dx_psi1 + dt_psi[1] * dx_psi0) / denominator

                #TODO: check if this factor of 2 is actually supposed to be there (Vidar 05.12.23)

                print(region_to_ignore.shape)
                print(Vx.shape)
                print(Vy.shape)

                Vx[region_to_ignore] = 0
                Vy[region_to_ignore] = 0

                return [Vx, Vy]

        elif self.dim == 3:
            if len(psi) == 2:
                # Parameters to exclude region
                threshold = 0.4

                psi_f = [sp.fft.fftn(psi[0]), sp.fft.fftn(psi[1])]

                dx_psi0 = sp.fft.ifftn(self.dif[0] * psi_f[0])
                dy_psi0 = sp.fft.ifftn(self.dif[1] * psi_f[0])
                dz_psi0 = sp.fft.ifftn(self.dif[2] * psi_f[0])

                dx_psi1 = sp.fft.ifftn(self.dif[0] * psi_f[1])
                dy_psi1 = sp.fft.ifftn(self.dif[1] * psi_f[1])
                dz_psi1 = sp.fft.ifftn(self.dif[2] * psi_f[1])


                denominator = 2 * (dx_psi1 ** 2 * (dy_psi0 ** 2 + dz_psi0 ** 2) + (
                            dy_psi1 * dz_psi0 - dy_psi0 * dz_psi1) ** 2 -
                                   2 * dx_psi0 * dx_psi1 * (dy_psi0 * dy_psi1 + dz_psi0 * dz_psi1) +
                                   dx_psi0 ** 2 * (dy_psi1 ** 2 + dz_psi1 ** 2))

                Vx = -2 * np.real((
                          (-dx_psi1 * (dy_psi0 * dy_psi1 + dz_psi0 * dz_psi1) +
                           dx_psi0 * (dy_psi1 ** 2 + dz_psi1 ** 2)) * dt_psi[0] +
                          (dx_psi1 * (dy_psi0 ** 2 + dz_psi0 ** 2) -
                           dx_psi0 * (dy_psi0 * dy_psi1 + dz_psi0 * dz_psi1)) * dt_psi[1]
                                  )/denominator)
                Vy = -2 * np.real((
                          (dx_psi1 ** 2 * dy_psi0 - dx_psi0 * dx_psi1 * dy_psi1 +
                           dz_psi1 * (-dy_psi1 * dz_psi0 + dy_psi0 * dz_psi1)) * dt_psi[0] +
                          (-dx_psi0 * dx_psi1 * dy_psi0 + dx_psi0 ** 2 * dy_psi1 +
                           dz_psi0 * (dy_psi1 * dz_psi0 - dy_psi0 * dz_psi1)) * dt_psi[1]
                                  ) / denominator)

                Vz = -2 * np.real((
                          ((dx_psi1 ** 2 + dy_psi1 ** 2) * dz_psi0 - (
                                      dx_psi0 * dx_psi1 + dy_psi0 * dy_psi1) * dz_psi1) * dt_psi[0] +
                          (-(dx_psi0 * dx_psi1 + dy_psi0 * dy_psi1) * dz_psi0 + (
                                      dx_psi0 ** 2 + dy_psi0 ** 2) * dz_psi1) * dt_psi[1]
                                ) / denominator)

                denominator_max = np.max(abs(denominator))
                region_to_ignore = abs(denominator) < threshold * denominator_max

                Vx[region_to_ignore] = 0
                Vy[region_to_ignore] = 0
                Vz[region_to_ignore] = 0

                return [Vx, Vy, Vz]

    def calc_defect_current_density(self, psi, dt_psi, psi_0=0):
        """
        Calculates the conserved current of the superfluid density
        Args:
            psi (numpy.ndarray) the vector field that we find the density of
            dt_psi (numpy.ndarray) the time derivative of psi
            psi_0 (floar or numpy.ndarray, optional) the equilibrium state
        :returns
            (array) components of the conserved current
        """
        if self.dim == 2:
            if len(psi) == 2:
                psi_f = [sp.fft.fftn(psi[0]), sp.fft.fftn(psi[1])]

                dx_psi0 = sp.fft.ifftn(self.dif[0] * psi_f[0])
                dy_psi1 = sp.fft.ifftn(self.dif[1] * psi_f[1])
                dx_psi1 = sp.fft.ifftn(self.dif[0] * psi_f[1])
                dy_psi0 = sp.fft.ifftn(self.dif[1] * psi_f[0])

                Jx = -  np.real(dt_psi[0] * dy_psi1 - dt_psi[1] * dy_psi0) / (psi_0 * np.pi)
                Jy = - np.real(-dt_psi[0] * dx_psi1 + dt_psi[1] * dx_psi0) / (psi_0 * np.pi)

                return [Jx, Jy]

    def calc_delta_function(self, psi, psi0=1):
        """
        Calculate the delta function for a given wavefunction.

        Parameters:
            psi (list): The wavefunction.
            psi0 (float): The width of the wavefunction. Default is 1.

        Returns:
            float: The value of the delta function.
        """
        width = psi0 / 2
        n = len(psi)
        if self.dim == 2:
            if n == 2:
                psi2 = psi[0] ** 2 + psi[1] ** 2
                return 1 / (2 * np.pi * width ** 2) * np.exp(-psi2 / (2 * width ** 2))

    def calc_region_disk(self, position, radius):

        if self.dim == 2:
            rx2m = (self.x - position[0] - self.xmax) ** 2
            rx2 = (self.x - position[0]) ** 2
            rx2p = (self.x - position[0] + self.xmax) ** 2
            rx2 = np.min(np.stack((rx2m, rx2, rx2p)), axis=0).reshape((self.xRes, 1))

            ry2m = (self.y - position[1] - self.ymax) ** 2
            ry2 = (self.y - position[1]) ** 2
            ry2p = (self.y - position[1] + self.ymax) ** 2
            ry2 = np.min(np.stack((ry2m, ry2, ry2p)), axis=0).reshape((1, self.yRes))

            return rx2 + ry2 <= radius ** 2

        else:
            raise Exception("Not valid for other dimensions.")

    def calc_region_ball(self, position, radius):
        """
        Calculates a boolean array indicating whether a point is within a ball of a given radius.
        """
        if self.dim == 3:
            # This code ensures that the region is periodic
            rx2m = (self.x - position[0] - self.xmax) ** 2
            rx2 = (self.x - position[0]) ** 2
            rx2p = (self.x - position[0] + self.xmax) ** 2
            rx2 = np.min(np.stack((rx2m, rx2, rx2p)), axis=0).reshape((self.xRes, 1, 1))

            ry2m = (self.y - position[1] - self.ymax) ** 2
            ry2 = (self.y - position[1]) ** 2
            ry2p = (self.y - position[1] + self.ymax) ** 2
            ry2 = np.min(np.stack((ry2m, ry2, ry2p)), axis=0).reshape((1, self.yRes, 1))

            rz2m = (self.z - position[2] - self.zmax) ** 2
            rz2 = (self.z - position[2]) ** 2
            rz2p = (self.z - position[2] + self.zmax) ** 2
            rz2 = np.min(np.stack((rz2m, rz2, rz2p)), axis=0).reshape((1, 1, self.zRes))

            return rx2 + ry2 + rz2 <= radius ** 2

        else:
            raise Exception("Not valid for other dimensions.")

    def calc_region_cylinder(self, position, radius, normal_vector, height):
        """
        Calculates a boolean array indicating whether a point is within a cylinder of a given radius and height.
        """
        if self.dim == 3:
            t = normal_vector / np.linalg.norm(np.array(normal_vector))

            rx = (self.x - position[0]).reshape((self.xRes, 1, 1))
            rx[rx > self.xmax / 2] = rx[rx > self.xmax / 2] - self.xmax
            rx[rx < -self.xmax / 2] = rx[rx < -self.xmax / 2] + self.xmax

            ry = (self.y - position[1]).reshape((1, self.yRes, 1))
            ry[ry > self.ymax / 2] = ry[ry > self.ymax / 2] - self.ymax
            ry[ry < -self.ymax / 2] = ry[ry < -self.ymax / 2] + self.ymax

            rz = (self.z - position[2]).reshape((1, 1, self.zRes))
            rz[rz > self.zmax / 2] = rz[rz > self.zmax / 2] - self.zmax
            rz[rz < -self.zmax / 2] = rz[rz < -self.zmax / 2] + self.zmax

            zt = rx * t[0] + ry * t[1] + rz * t[2]

            # Project to perpendicular plane vector
            Rt2 = (rx - zt * t[0]) ** 2 + (ry - zt * t[1]) ** 2 + (rz - zt * t[2]) ** 2

            return (zt ** 2 <= height ** 2) & (Rt2 <= radius ** 2)



        else:
            raise Exception("Not valid for other dimensions.")

    def calc_integrate_field(self, field, region=None):
        """
        Calculates the integrated field value within a specified region.

        Parameters:
            field (numpy.ndarray): The field array.
            index (tuple, optional): The indices of the center point in the field. Defaults to None.
            radius (float, optional): The radius of the region. Defaults to None.

        Returns:
            tuple or float: If index is provided, returns a tuple containing the integrated field value
                            within the region and a boolean array indicating the region. If index is None,
                            returns the integrated field value within the entire field.

        Raises:
            Exception: If the dimension of the field is not 2.
        """

        if region is None:
            return np.sum(field) * self.dV
        else:
            return np.sum(field[region]) * self.dV

    def calc_integrating_factors_f_and_solver(self, omega_f, method):

        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop

        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise Exception("This method is not implemented.")

        return integrating_factors_f, solver

    def calc_evolution_integrating_factors_ETD2RK(self, omega_f, tol=10 ** (-4)):
        """
        Calculates integrating factors for ETD2RK
         parameters
            omega_f (numpy.ndarray): the value of omega_f
            tol (float, optional): tolerance for when to expand the integrating factors that divide by omega
        :return:
            list: the list of integrating factors
        """
        integrating_factors_f = [0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = (If1 - 1) / omega_f
        integrating_factors_f[1][np.abs(omega_f) < tol] = self.dt

        integrating_factors_f[2] = 1 / (self.dt * omega_f ** 2) * (If1 - 1 - omega_f * self.dt)
        integrating_factors_f[2][np.abs(omega_f) < tol] = self.dt / 2
        return integrating_factors_f

    def calc_evolution_integrating_factors_ETD4RK(self, omega_f, tol=10 ** (-4)):
        """
        Calculate the evolution integrating factors using the ETDRK4 method.

        Parameters:
            omega_f (numpy.ndarray): The value of omega_f.
            tol (float,optional): tolerance for when to expand the integrating factors that divide by omega
        Returns:
            list: The list of integrating factors.
        """
        integrating_factors_f = [0, 0, 0, 0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt / 2)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = (If1 - 1) / omega_f
        integrating_factors_f[1][np.abs(omega_f) < tol] = self.dt / 2

        integrating_factors_f[2] = np.exp(omega_f * self.dt)

        integrating_factors_f[3] = 1 / (self.dt ** 2 * omega_f ** 3) \
                                   * (-4 - self.dt * omega_f + If1 ** 2 * (
                4 - 3 * self.dt * omega_f + self.dt ** 2 * omega_f ** 2))

        integrating_factors_f[3][np.abs(omega_f) < tol] = self.dt / 6

        integrating_factors_f[4] = 1 / (self.dt ** 2 * omega_f ** 3) \
                                   * (2 + self.dt * omega_f + If1 ** 2 * (-2 + self.dt * omega_f))

        integrating_factors_f[4][np.abs(omega_f) < tol] = self.dt / 3

        integrating_factors_f[5] = 1 / (self.dt ** 2 * omega_f ** 3) \
                                   * (-4 - 3 * self.dt * omega_f - self.dt ** 2 * omega_f ** 2 + If1 ** 2 * (
                4 - self.dt * omega_f))

        integrating_factors_f[5][np.abs(omega_f) < tol] = self.dt / 6

        return integrating_factors_f

    def calc_defect_nodes(self, defect_density,
                          charge_tolerance=None,
                          integration_radius=None):
        """
        Calculate the positions and charges of defect nodes based on the defect density.
        Input:
            defect_density (numpy.ndarray): The defect density field. A positive scalar field to be integrated.
        Returns:
            list: A list of dictionaries representing the defect nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the defect node in the defect density array.
                  - 'position': The position of the defect node
        """

        if not (defect_density>=0).all():
            raise Exception("Defect density must be a positive real scalar field.")

        defect_nodes = []

        if self.dim == 2:
            if charge_tolerance is None:
                charge_tolerance = 0.2
            if integration_radius is None:
                integration_radius = self.a0
            
            #Auxiliary functions
            def calc_region(defect_density_max_index,radius):
                return self.calc_region_disk(position=[
                    self.x.flatten()[defect_density_max_index[0]],
                    self.y.flatten()[defect_density_max_index[1]]], 
                    radius=radius)

            def calc_position_from_region(defect_density,region_to_integrate):
                x = np.sum(region_to_integrate * defect_density * self.x) / np.sum(region_to_integrate * defect_density)
                y = np.sum(region_to_integrate * defect_density * self.y) / np.sum(region_to_integrate * defect_density)
                return [x,y]

        elif self.dim == 3:
            if charge_tolerance is None:
                charge_tolerance = 0.5*self.a0**2
                # print("charge tolerance",charge_tolerance)
            if integration_radius is None:
                integration_radius = 2*self.a0

            #Auxiliary functions
            def calc_region(defect_density_max_index,radius):
                return self.calc_region_ball(position=[
                        self.x.flatten()[defect_density_max_index[0]],
                        self.y.flatten()[defect_density_max_index[1]],
                        self.z.flatten()[defect_density_max_index[2]]], 
                        radius=radius)
            
            def calc_position_from_region(defect_density,region_to_integrate):
                x = np.sum(region_to_integrate * defect_density * self.x) / np.sum(region_to_integrate * defect_density)
                y = np.sum(region_to_integrate * defect_density * self.y) / np.sum(region_to_integrate * defect_density)
                z = np.sum(region_to_integrate * defect_density * self.z) / np.sum(region_to_integrate * defect_density)
                return [x,y,z]

        #Region to search for defect nodes
        region_to_search = np.ones(self.dims)

        # Calculate the point where defect density is largest
        defect_density_max_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

        # Integrate the defect density around this point (i.e. in a disk/ball around)
        region_to_integrate = calc_region(defect_density_max_index,
                                        radius=integration_radius)
        
        region_to_exclude_from_search = calc_region(defect_density_max_index,
                                        radius=2*integration_radius)
        
        region_to_search[region_to_exclude_from_search] = 0
        
        charge = self.calc_integrate_field(defect_density, region_to_integrate)

        while charge > charge_tolerance:
            defect_node = {}
            defect_node['position_index'] = defect_density_max_index

            defect_node['position'] = calc_position_from_region(defect_density,region_to_integrate)

            defect_nodes.append(defect_node)

            # print("charge:", charge)
            # self.plot_field(defect_density)
            # plt.show()

            defect_density[region_to_integrate] = 0

            defect_density_max_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

            region_to_integrate = calc_region(defect_density_max_index,
                                        radius=integration_radius)
            
            region_to_exclude_from_search = calc_region(defect_density_max_index,
                                        radius=2*integration_radius)
            
            region_to_search[region_to_exclude_from_search] = 0
            
            charge = self.calc_integrate_field(defect_density, region_to_integrate)

        return defect_nodes

    ## Time evolution function
    def evolve_ETD2RK_loop(self, integrating_factors_f, non_linear_evolution_function_f, field, field_f):
        """
        Evolves the given field using the ETD2RK scheme with a loop.

        Parameters:
            integrating_factors_f (list): A list of three integrating factors.
            non_linear_evolution_function_f (function): A function that calculates the non-linear evolution of the field
                    and returns the fourier transform.
            field (ndarray): The initial field to be evolved.
            field_f (ndarray): The Fourier transform of the initial field.

        Returns:
            tuple: A tuple containing the evolved field and the predicted field in Fourier space.
        """
        """
        Don't think this is needed. self.time is now initialized in base and all timedependence should be through this
        if t==None:
            def N_f(field,t):
                return non_linear_evolution_function_f(field)
            t=0

        else:
            def N_f(field,t):
                return non_linear_evolution_function_f(field,t)
        """
        t_0 = self.time
        N0_f = non_linear_evolution_function_f(field)

        a_f = integrating_factors_f[0] * field_f + integrating_factors_f[1] * N0_f
        a = sp.fft.ifftn(a_f, axes=(range(-self.dim, 0)))
        self.time = t_0 + self.dt
        N_a_f = non_linear_evolution_function_f(a)
        field_f = a_f + integrating_factors_f[2] * (N_a_f - N0_f)
        field = sp.fft.ifftn(field_f, axes=(range(-self.dim, 0)))

        return field, field_f

    def evolve_ETD4RK_loop(self, integrating_factors_f, non_linear_evolution_function_f, field, field_f):

        """
         Evolves the given field using the ETD4RK scheme with a loop.

         Parameters:
             integrating_factors_f (list): A list of five integrating factors.
             non_linear_evolution_function_f (function): A function that calculates the non-linear evolution of the field.
             field (ndarray): The initial field to be evolved.
             field_f (ndarray): The Fourier transform of the initial field.

         Returns:
             tuple: A tuple containing the evolved field and the predicted field in Fourier space.
         """
        N_0f = non_linear_evolution_function_f(field)

        t_0 = self.time
        self.time = t_0 + self.dt / 2

        a_f = field_f * integrating_factors_f[0] + N_0f * integrating_factors_f[1]
        a = sp.fft.ifftn(a_f, axes=(range(-self.dim, 0)))
        N_a = non_linear_evolution_function_f(a)

        b_f = field_f * integrating_factors_f[0] + N_a * integrating_factors_f[1]
        b = sp.fft.ifftn(b_f, axes=(range(-self.dim, 0)))
        N_b = non_linear_evolution_function_f(b)

        self.time = t_0 + self.dt
        c_f = a_f * integrating_factors_f[0] + (2 * N_b - N_0f) * integrating_factors_f[1]
        c = sp.fft.ifftn(c_f, axes=(range(-self.dim, 0)))
        N_c = non_linear_evolution_function_f(c)

        field_f = field_f * integrating_factors_f[2] + N_0f * integrating_factors_f[3] \
                  + 2 * (N_a + N_b) * integrating_factors_f[4] + N_c * integrating_factors_f[5]

        field = sp.fft.ifftn(field_f, axes=(range(-self.dim, 0)))

        return field, field_f

    def evolve_ETD2RK_loop_test(self, integrating_factors_f, non_linear_evolution_function_f, field, field_f):
        # TODO this can eventualy be removed
        N_0 = non_linear_evolution_function_f(field)
        field_f = field_f * integrating_factors_f[0] + N_0 * integrating_factors_f[1]
        temp = sp.fft.ifftn(field_f, axes=(range(-self.dim, 0)))
        N_1 = non_linear_evolution_function_f(temp) - N_0
        field_f += N_1 * integrating_factors_f[2]
        field = sp.fft.ifftn(field_f, axes=(range(-self.dim, 0)))
        return field, field_f

    # PLOTTING FUNCTIONS

    def plot_angle_field(self, field, ax=None, colorbar=True):
        """
        Plot the angle field.
        Parameters:
            field (array-like): The angle field values.
            ax (matplotlib.axes.Axes, optional): The axes to plot the angle field on. If not provided, a new subplot will be created.
        Returns:
            None
        """

        if self.dim == 2:

            if ax is None:
                ax = plt.gca()

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            custom_colormap = tool_colormap_angle()

            mesh = ax.pcolormesh(X, Y, field, shading='auto', cmap=custom_colormap, vmin=-np.pi, vmax=np.pi)
            if colorbar:
                cbar = plt.colorbar(mesh)  # To add a colorbar on the side
                cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
                cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])
            # ax.title("Angle field")
            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_aspect('equal')

            return ax

        elif self.dim == 3:

            if ax == None:
                plt.figure()
                ax = plt.gcf().add_subplot(111, projection='3d')

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            cmap = tool_colormap_angle()

            field_min = np.min(field)
            field_max = np.max(field)

            for angle in [-2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3]:

                if field_min < angle < field_max:
                    field_to_plot = field.copy()
                    field_to_plot[field < angle - 1] = float('nan')
                    field_to_plot[field > angle + 1] = float('nan')

                    verts, faces, _, _ = marching_cubes(field_to_plot, angle)

                    ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                                    color=cmap((angle + np.pi) / (2 * np.pi)))

            field = np.mod(field, 2 * np.pi)

            field_to_plot = field.copy()
            field_to_plot[field < np.pi - 1] = float('nan')
            field_to_plot[field > np.pi + 1] = float('nan')

            verts, faces, _, _ = marching_cubes(field_to_plot, np.pi)

            ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                            color=cmap(0))

            ax.set_xlim3d(self.xmin, self.xmax-self.dx)
            ax.set_ylim3d(self.ymin, self.ymax-self.dy)
            ax.set_zlim3d(self.zmin, self.zmax-self.dz)
            ax.set_aspect('equal')

    def plot_field(self, field, ax=None, colorbar=True, colormap=None, cmax=None, cmin=None,
                   clims = None,
                   number_of_layers=1, hold=False, cmap_symmetric=None, layer_values=None):
        """
        Plots the given field.
        Parameters:
            field (array-like): The field to be plotted.
            ax (Axes, optional): The axes object to plot on. If None, a new figure and axes will be created.
            colorbar (bool, optional): Whether to include a colorbar in the plot. Default is True.
            colormap (Colormap, optional): The colormap to use for the plot. If None, a default colormap will be used.
            cmax (float, optional): The maximum value for the colorbar. If None, the maximum value of the field will be used.
            cmin (float, optional): The minimum value for the colorbar. If None, the minimum value of the field will be used.
            number_of_layers (int, optional): The number of layers to plot for a 3D field. Default is 1.
            hold (bool, optional): Whether to clear the axes before plotting. Default is False.
            cmap_symmetric (bool, optional): Whether to make the colormap symmetric. Default is True.
        Returns:
            ax (Axes): The axes object containing the plot.
        """



        if self.dim == 1:
            if ax == None:
                plt.clf()
                ax = plt.gca()

            ax.plot(self.x/self.a0, field)
            ax.set_xlabel('$x/a_0$')
            ax.grid(True)


        if field.dtype == bool:
            field = field.astype(float)

        if self.dim == 2:

            if ax == None:
                plt.clf()
                ax = plt.gca()

            if colormap is None:
                cmap = plt.get_cmap('viridis')

            elif colormap == 'bluewhitered':
                cmap = tool_colormap_bluewhitered()
                # Set symmetric if not specified otherwise
                if cmap_symmetric is None:
                    print("Setting the colormap to be symmetric with use of bluewhitered since not otherwiese specified.")
                    cmap_symmetric = True

            elif colormap == 'parula':
                cmap = tool_colormap_parula()

            else:
                cmap = plt.get_cmap(colormap)

            if cmap_symmetric is None:
                cmap_symmetric = False

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            pcm = ax.pcolormesh(X / self.a0, Y / self.a0, field, shading='gouraud', cmap=cmap)
            ax.set_aspect('equal')

            if clims is not None:
                cmin = clims[0]
                cmax = clims[1]
            if cmin is not None:
                pcm.set_clim(vmin=cmin)
            if cmax is not None:
                pcm.set_clim(vmax=cmax)

            if cmap_symmetric:
                cmax = abs(field).max()
                cmin = -cmax
                pcm.set_clim(vmin=cmin, vmax=cmax)

            if colorbar:
                cbar = plt.colorbar(pcm, ax=ax)

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_aspect('equal')

            return ax

        elif self.dim == 3:

            if ax == None:
                plt.clf()
                ax = plt.gcf().add_subplot(111, projection='3d')

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            field_min = np.min(field)
            field_max = np.max(field)

            if clims is None:
                cmin = field_min
                cmax = field_max
            else:
                cmin = clims[0]
                cmax = clims[1]

            if layer_values is None:
                layer_values = np.linspace(cmin, cmax, number_of_layers + 2)
            else:
                layer_values = np.concatenate([[-np.inf], layer_values, [np.inf]])

            # print(layer_values)
            if colormap is None:
                cmap = plt.get_cmap('viridis')
            else:
                cmap = plt.get_cmap(colormap)

            if field_min < layer_values[1] < field_max:
                verts, faces, _, _ = marching_cubes(field, layer_values[1])
                ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                            color=cmap(layer_values[1] / field_max))

            for layer_value in layer_values[2:-1]:

                if field_min < layer_value < field_max:
                    verts, faces, _, _ = marching_cubes(field, layer_value)
                    ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                                color=cmap(layer_value / field_max))

            ax.set_aspect('equal')
            if colorbar:
                sm = plt.cm.ScalarMappable(cmap=cmap)
                sm.set_clim(cmin, cmax)
                plt.colorbar(sm, ax=ax)

            ax.set_xlim3d(self.xmin, self.xmax-self.dx)
            ax.set_ylim3d(self.ymin, self.ymax-self.dy)
            ax.set_zlim3d(self.zmin, self.zmax-self.dz)
            ax.set_aspect('equal')
            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')


            return ax

    def plot_fourier_field(self, field_f, ax=None):
        """
            Plot a Fourier field.

            Parameters:
                field_f (ndarray): The Fourier field to be plotted.
                ax (Axes3D, optional): The matplotlib 3D axis to be used for plotting. If not provided, a new axis will be created.

            Returns:
                None
            """
        field_f = np.fft.fftshift(field_f)

        if ax == None:
            ax = plt.gcf().add_subplot(111, projection='3d')

        if self.dim == 2:
            rho = np.abs(field_f)
            theta = np.angle(field_f)

            Kx, Ky = np.meshgrid(self.k[0], self.k[1], indexing='ij')

            Kx = np.fft.fftshift(Kx)
            Ky = np.fft.fftshift(Ky)

            custom_colormap = tool_colormap_angle()

            # Get the colors from a colormap (e.g., hsv, but you can choose any other)
            colors = plt.cm.hsv((theta + np.pi) / (2 * np.pi))  # Normalizing theta to [0, 1]
            surf = ax.plot_surface(Kx, Ky, rho, facecolors=colors, shade=True)

            # mappable = plt.cm.ScalarMappable(cmap=custom_colormap)
            # mappable.set_array([])
            # mappable.set_clim(-np.pi, np.pi)
            # cbar = plt.colorbar(mappable, ax=ax)
            # cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
            # cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

            # plt.title("Angle field")
            # plt.xlabel("X-axis")
            # plt.ylabel("Y-axis")

    def plot_complex_field(self, complex_field, ax=None, plot_method=None, colorbar=False):
        """
        Plot a complex field.

        Parameters:
            complex_field (numpy.ndarray): The complex field to plot.
            ax (matplotlib.axes.Axes, optional): The matplotlib axes on which to plot the field.
                If not provided, a new 3D axes will be created.

        Raises:
            Exception: If the dimension of the field is not 2.

        Returns:
            None
        """



        if self.dim == 2:
            if plot_method is None:
                plot_method = 'phase_angle'

            if plot_method == '3Dsurface':
                if ax == None:
                    ax = plt.gcf().add_subplot(111, projection='3d')

                X, Y = np.meshgrid(self.x, self.y, indexing='ij')

                custom_colormap = tool_colormap_angle()

                # Get the colors from a colormap (e.g., hsv, but you can choose any other)
                colors = plt.cm.hsv((theta + np.pi) / (2 * np.pi))  # Normalizing theta to [0, 1]

                surf = ax.plot_surface(X, Y, rho, facecolors=colors)
            elif plot_method == 'phase_angle':
                if ax == None:
                    plt.clf()
                    ax = plt.gca()

                X, Y = np.meshgrid(self.x, self.y, indexing='ij')

                rho = np.abs(complex_field)
                theta = np.angle(complex_field)

                rho_normalized = rho / np.max(rho)
                custom_colormap_phase = tool_colormap_angle()


                # Create a new colormap for magnitudeW
                # Starting from white (for zero magnitude) to the full color of the phase
                # custom_colormap_mag = mcolors.LinearSegmentedColormap.from_list(
                #     'MagnitudeColorMap',
                #     [(1, 1, 1, 0), custom_colormap_phase(1.0)],
                #     N=256
                # )

                # Calculate colors based on magnitude and phase
                #colors = custom_colormap_phase(theta)
                #colors[..., 3] = rho_normalized  # Set the alpha channel according to the magnitude

                mesh = ax.pcolormesh(X, Y, theta, shading='auto', cmap=custom_colormap_phase, vmin=-np.pi, vmax=np.pi)
                mesh.set_alpha(rho_normalized)
                
                #mesh.set_array(None)  # Avoids warning
                #mesh.set_edgecolor('face')
                #mesh.set_facecolor(colors)  # Use the calculated colors

            if colorbar:
                mappable = plt.cm.ScalarMappable(cmap=custom_colormap)
                mappable.set_array([])
                mappable.set_clim(-np.pi, np.pi)
                cbar = plt.colorbar(mappable, ax=ax)
                cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
                cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

            ax.set_xlabel("$x/a_0$")
            ax.set_ylabel("$y/a_0$")
            ax.set_aspect('equal')
            
            return ax

        elif self.dim == 3:

            if plot_method is None:
                plot_method = 'phase_blob'

            if ax == None:
                plt.clf()
                ax = plt.gcf().add_subplot(111, projection='3d')

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            rho = np.abs(complex_field)
            rho_normalized = rho / np.max(rho)
            theta = np.angle(complex_field)

            cmap = tool_colormap_angle()

            if plot_method == 'phase_angle':
                

                for angle in [-2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3]:
                    field_to_plot = theta.copy()
                    field_to_plot[theta < angle - 1] = float('nan')
                    field_to_plot[theta > angle + 1] = float('nan')
                    field_to_plot[rho_normalized < 0.01] = float('nan')

                    if np.nanmin(field_to_plot) < angle < np.nanmax(field_to_plot):

                        verts, faces, _, _ = marching_cubes(field_to_plot, angle)

                        ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                                        color=cmap((angle + np.pi) / (2 * np.pi)))

                theta = np.mod(theta, 2 * np.pi)

                field_to_plot = theta.copy()
                field_to_plot[theta < np.pi - 1] = float('nan')
                field_to_plot[theta > np.pi + 1] = float('nan')
                field_to_plot[rho_normalized < 0.01] = float('nan')

                if np.nanmin(field_to_plot) < np.pi < np.nanmax(field_to_plot):

                    verts, faces, _, _ = marching_cubes(field_to_plot, np.pi)

                    ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, self.ymin+verts[:, 1]*self.dy, faces, self.zmin+verts[:, 2]*self.dz, alpha=0.5,
                                color=cmap(0))
            
            elif plot_method == 'phase_blob':

                if np.nanmin(rho_normalized)<0.5<np.nanmax(rho_normalized):
                    verts, faces, _, _ = marching_cubes(rho_normalized, 0.5)

                    # Calculate the centroids of each triangle
                    centroids = np.mean(verts[faces], axis=1)

                    # Assuming theta is defined on the same grid as rho
                    x, y, z = np.mgrid[0:rho_normalized.shape[0], 0:rho_normalized.shape[1], 0:rho_normalized.shape[2]]

                    # Flatten the grid for interpolation
                    points = np.c_[x.ravel(), y.ravel(), z.ravel()]
                    theta_values = theta.ravel()

                    # Interpolate theta at the vertices positions
                    theta_verts = sp.interpolate.griddata(points, theta_values, centroids, method='nearest')

                    # Normalize theta values for color mapping
                    theta_normalized = (theta_verts + np.pi) / (2*np.pi)

                    # Map normalized theta values to colors
                    colors = cmap(theta_normalized)

                    # print("Colors shape:", colors.shape)
                    # print(colors)

                    ax.plot_trisurf(self.xmin+verts[:, 0]*self.dx, 
                                    self.ymin+verts[:, 1]*self.dy, 
                                    faces, 
                                    self.zmin+verts[:, 2]*self.dz, 
                                    facecolor=colors, antialiased=False)

            ax.set_xlim3d(self.xmin, self.xmax-self.dx)
            ax.set_ylim3d(self.ymin, self.ymax-self.dy)
            ax.set_zlim3d(self.zmin, self.zmax-self.dz)

            ax.set_aspect('equal')

            ax.set_xlabel("$x/a_0$")
            ax.set_ylabel("$y/a_0$")
            ax.set_zlabel("$z/a_0$")
            ax.set_aspect('equal')


        else:
            raise Exception("This plotting function not yet configured for other dimension")

    def plot_field_in_plane(self, field, normal_vector=[0,1,0], position=None, ax=None,
                            colorbar=True, colormap = 'winter'):

        if self.dim != 3:
            raise Exception("This plotting function not yet configured for other dimensions")

        if position is None:
            position = self.rmid

        if ax is None:
            plt.clf()
            ax = plt.gcf().add_subplot(111, projection='3d')

        unit_vector1, unit_vector2, t = tool_create_orthonormal_triad(normal_vector)

        sizes = [self.xmax-self.xmin-self.dx,
                 self.ymax-self.ymin-self.dy,
                 self.zmax-self.zmin-self.dz]
        max_size = max(sizes)
        dr = [self.dx, self.dy, self.dz]
        min_step = min(dr)
        linear_mesh = np.arange(-1*max_size,1*max_size, min_step)

        #Mesh containing x- and y-coordinates on the plane
        [X0,Y0] = np.meshgrid(linear_mesh, linear_mesh, indexing='ij')

        #Mesh containing the true x-, y- and z-coordinates on the plane
        R = [position[i] + unit_vector1[i] * X0 + unit_vector2[i] * Y0 for i in range(3)]

        #Remove excess points from R outside of the domain
        points_to_exclude = np.zeros_like(X0)
        points_to_exclude += (R[0] > (self.xmax-self.dx)) + (R[0] < self.xmin)
        points_to_exclude += (R[1] > (self.ymax-self.dy)) + (R[1] < self.ymin)
        points_to_exclude += (R[2] > (self.zmax-self.dz)) + (R[2] < self.zmin)
        points_to_exclude = points_to_exclude.astype(bool)

        # print("xmin:", self.xmin)
        # print("Lmax:", self.xmax-self.dx)
        # print("ymin:", self.ymin)
        # print("Lmax:", self.ymax-self.dy)
        # print("zmin:", self.zmin)
        # print("Lmax:", self.zmax-self.dz)
        for i in range(3):
            R[i][points_to_exclude] = 0

        #For each point in R, find:
        #   1. The true index (float) of the point in the matrix closest to the point in R
        rmin = [self.xmin, self.ymin, self.zmin]
        Ri = [(R[i]-rmin[i])/dr[i] for i in range(3)]

        # 2. The lower and upper limit of the integer index of the point in the matrix closest to the point in R
        RiLH = []
        Res = [self.xRes,self.yRes,self.zRes]
        for i in range(3):
            Ri_floor = np.floor(Ri[i])

            Ri_ceil = np.ceil(Ri[i])

            # Adjusting floor and ceil for integer values of Ri
            is_integer = Ri[i] == Ri[i].astype(int)
            # print(is_integer)
            # Decrease floor by 1 if Ri is not at lower boundary
            Ri_floor[is_integer & (Ri[i] > 0)] -= 1

            # Increase ceil by 1 if Ri is not at upper boundary
            Ri_ceil[is_integer & (Ri[i] == 0)] += 1

            RiLH.append([Ri_floor, Ri_ceil])


        # 3. The index distances between the point in R and the point in the matrix closest to it
        dRi = [[Ri[i] - RiLH[i][0], RiLH[i][1] - Ri[i]] for i in range(3)]

        field_on_plane = np.zeros_like(X0)
        resolution_on_plane = len(linear_mesh)

        for ix in [0, 1]:
            for iy in [0, 1]:
                for iz in [0, 1]:
                    field_on_plane= field_on_plane \
                                           + field[RiLH[0][ix].astype(int),
                                                    RiLH[1][iy].astype(int),
                                                    RiLH[2][iz].astype(int)] * \
                                               (1 - dRi[0][ix]) * \
                                               (1 - dRi[1][iy]) * \
                                               (1 - dRi[2][iz])

        # field_on_plane[points_to_exclude] = np.nan
        for i in range(3):
            R[i][points_to_exclude] = np.nan

        cmin = np.nanmin(field_on_plane)
        cmax = np.nanmax(field_on_plane)

        # Normalize the field data
        norm = mcolors.Normalize(vmin=cmin, vmax=cmax, clip=True)

        if colormap == 'angle':
            cmap = tool_colormap_angle()
        else:
            cmap = plt.get_cmap(colormap)

        # Apply the colormap to your field datav
        field_on_plane_colors = cmap(norm(field_on_plane))

        # Set NaNs to transparent
        field_on_plane_colors[np.isnan(field_on_plane)] = [0, 0, 0, 0]

        ax.plot_surface(R[0] / self.a0, R[1] / self.a0, R[2] / self.a0, facecolors=field_on_plane_colors,
                        rstride=1,cstride=1)


        ax.set_xlim3d(self.xmin/self.a0, (self.xmax - self.dx)/self.a0)
        ax.set_ylim3d(self.ymin/self.a0, (self.ymax - self.dy)/self.a0)
        ax.set_zlim3d(self.zmin/self.a0, (self.zmax - self.dz)/self.a0)
        ax.set_aspect('equal')
        ax.set_xlabel('$x/a_0$')
        ax.set_ylabel('$y/a_0$')
        ax.set_zlabel('$z/a_0$')

    def plot_angle_field_in_plane(self, angle_field, colorbar = True):

        self.plot_field_in_plane(angle_field, colorbar=False)


        if colorbar:
            sm = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
            sm.set_clim(-np.pi, np.pi)
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_ticks(np.array([0,1/6,2/6,3/6,4/6,5/6,1]))
            cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

    def plot_vector_field(self, vector_field, ax=None, step=None):
        """
            Plots a vector field on a 2D grid.

            Parameters:
            vector_field (tuple): Tuple containing the x and y components of the vector field.
            ax (matplotlib.axes.Axes, optional): The axes on which to plot the vector field. If not provided, a new subplot will be created.

            Returns:
            None
            """

        if self.dim == 2:

            if ax == None:
                ax = plt.gcf().add_subplot(111)

            if step == None:
                step = 5

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            X_plot = X[::step, ::step]
            Y_plot = Y[::step, ::step]
            U_plot = vector_field[0][::step, ::step]
            V_plot = vector_field[1][::step, ::step]

            max_vector = np.max(np.sqrt(U_plot ** 2 + V_plot ** 2))
            print(max_vector)

            ax.quiver(X_plot, Y_plot, U_plot, V_plot, scale=25 * max_vector / step)

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_aspect('equal')
            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])

        elif self.dim == 3:

            if ax == None:
                ax = plt.gcf().add_subplot(111, projection='3d')

            if step is None:
                step = 2

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            X_plot = X[::step, ::step, ::step]
            Y_plot = Y[::step, ::step, ::step]
            Z_plot = Z[::step, ::step, ::step]
            U_plot = vector_field[0][::step, ::step, ::step]
            V_plot = vector_field[1][::step, ::step, ::step]
            W_plot = vector_field[2][::step, ::step, ::step]

            max_vector = np.max(np.sqrt(U_plot ** 2 + V_plot ** 2 + W_plot ** 2))

            U_plot = U_plot / max_vector
            V_plot = V_plot / max_vector
            W_plot = W_plot / max_vector

            ax.quiver(X_plot, Y_plot, Z_plot, U_plot, V_plot, W_plot, arrow_length_ratio=0.6)

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')
            ax.set_aspect('equal')
            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])
            ax.set_zlim([0, self.zmax-self.dz])

    def plot_field_velocity_and_director(self, field, velocity, director, ax=None, colorbar=True, colormap='viridis',
                                         cmax=None, cmin=None,
                                         number_of_layers=1, hold=False):
        """
        Plot the field, velocity, and director in the given axes.

        Parameters:
            field (ndarray): The field to be plotted.
            velocity (ndarray): The velocity to be plotted.
            director (ndarray): The director to be plotted.
            ax (Axes, optional): The axes to plot the field, velocity, and director on. If not provided, a new subplot will be created.
            colorbar (bool, optional): Whether to show the colorbar. Default is True.
            colormap (str, optional): The colormap to use for plotting the field. Default is 'viridis'.
            cmax (float, optional): The maximum value for the colorbar. If not provided, the maximum value of the field will be used.
            cmin (float, optional): The minimum value for the colorbar. If not provided, the minimum value of the field will be used.
            number_of_layers (int, optional): The number of layers in the plot. Default is 1.
            hold (bool, optional): Whether to hold the plot. Default is False.

        Returns:
            ax (Axes): The axes with the plotted field, velocity, and director.

        Raises:
            Exception: If the plotting function is not yet configured for dimensions other than 2.

        Note: streamplot assumes xy indexing and not ij. I think it is suficient
        just transpose the matrices before putting them in
        """
        if self.dim == 2:
            if ax == None:
                ax = plt.gcf().add_subplot(111)
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            cbar = ax.pcolormesh(X, Y, field, shading='gouraud', cmap=colormap)
            plt.colorbar(cbar)
            ax.streamplot(X.T, Y.T, (velocity[0]).T, (velocity[1]).T, color='w')
            ax.quiver(X, Y, director[0], director[1], headwidth=0, scale=50)
            ax.set_aspect('equal')
            return ax

        else:
            raise Exception("This plotting function not yet configured for other dimension")
