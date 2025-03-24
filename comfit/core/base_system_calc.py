from typing import Literal, Optional, Union

import numpy as np
import scipy as sp

from comfit.tool.tool_math_functions import tool_multinom

class BaseSystemCalc:
    """Calculation method for the base system class"""
    # CALCULATION FUNCTIONS

    # Calculation of angle fields for vortices of different types
    def calc_angle_field_single_vortex(
        self,
        position: Optional[list] = None,
        charge: int = 1
    ) -> np.ndarray:
        """Calculate the angle field due to a single vortex.
        
        Parameters
        ----------
        position : list, optional
            The position of the vortex. Defaults to None.
        charge : int, optional
            The charge of the vortex. Defaults to 1.
        
        Returns
        -------
        numpy.ndarray
            The angle field calculated for the vortex.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 2.
        """
        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]
        
        theta = charge * np.arctan2(self.y - position[1], self.x - position[0])

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_angle_field_vortex_dipole(
        self,
        dipole_vector: Optional[Union[list[float], tuple[float]]] = None,
        dipole_position: Optional[list[float]] = None
    ) -> np.ndarray:
        """Calculates the angle field for a double vortex system.

        Parameters
        ----------
        dipole_vector : list or tuple, optional
            The dipole vector. If None, defaults to [size_x/3, 0].
        dipole_position : list, optional
            The position of the center of mass of the dipole. If None, defaults to the middle of the grid.
            
        Returns
        -------
        numpy.ndarray
            The calculated angle field for the double vortex system.

        Raises
        ------
        Exception
            If the dimension of the system is not 2.
        """

        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if dipole_vector is None:
            dipole_vector = [self.size_x / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid
        
        theta = np.zeros(self.dims)
        theta += self.calc_angle_field_single_vortex(dipole_position - np.array(dipole_vector) / 2,
                                                            charge=-1)
        theta += self.calc_angle_field_single_vortex(dipole_position + np.array(dipole_vector) / 2, 
                                                            charge=1)

        # Convert the field to a complex field to make it fit the periodic boundary conditions
        amp = np.exp(1j * theta)

        # Filter the angle field
        width = 0.2 * np.min([self.size_x, self.size_y])
        radius = 0.4 * np.min([self.size_x, self.size_y])

        r2 = self.calc_distance_squared_to_point(self.rmid)
        filter = (1 + np.tanh((radius ** 2 - r2) / width ** 2)) / 2
        amp = amp * filter + (1 - filter)

        theta = np.angle(amp)

        # Roll the field so that the dipole center is at the desired position.
        Rx = round((self.rmid[0] - dipole_position[0]) / self.dx)
        theta = np.roll(theta, -Rx, axis=0)
        Ry = round((self.rmid[1] - dipole_position[1]) / self.dy)
        theta = np.roll(theta, -Ry, axis=1)

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_angle_field_vortex_ring(
            self,
            position: Optional[list[float]] = None,
            radius: Optional[float] = None,
            normal_vector: np.ndarray = [0, 0, 1]
            ) -> np.ndarray:
        """Calculate the angle field for a vortex ring.
        
        Parameters
        ----------
        position : list, optional
            The position of the vortex ring. If None, defaults to the middle of the grid.
        radius : float, optional
            The radius of the vortex ring. If None, defaults to 1/3 of the minimum system size.
        normal_vector : numpy.ndarray, optional
            The normal vector of the vortex ring. Defaults to [0, 0, 1].
            
        Returns
        -------
        numpy.ndarray
            The calculated angle field for the vortex ring.
        """
        if position is None:
            position = self.rmid

        if radius is None:
            radius = np.min([self.size_x, self.size_y, self.size_z]) / 3

        if radius > np.min([self.size_x, self.size_y, self.size_z]) / 3:
            print("Warning: The radius of the suggested vortex ring is large."
                  "This can cause unwanted boundary effects.")

        n = normal_vector / np.linalg.norm(np.array(normal_vector))
        [X, Y, Z] = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Add the vortices to the theta-field
        theta = 0
        position = np.array(position)


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
        width = 0.2 * np.min([self.size_x, self.size_y, self.size_z])
        radius = 0.4 * np.min([self.size_x, self.size_y, self.size_z])
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

    def calc_wavenums(
            self, 
            x: np.ndarray
            ) -> np.ndarray:
        """Calculates the wavenumbers corresponding to the input position vectors given by x.

        Parameters
        ----------
        x : numpy.ndarray
            1D array of x-positions

        Returns
        -------
        numpy.ndarray
            1D array of wavenumbers with all the modes for the given x-array,
            assuming periodicity from x[0] to x[0] over n intervals.

        Examples
        --------
        >>> x = np.array([-10, -5, 0, 5, 10])
        >>> k = instance_of_BaseSystem.calc_wavenums(self, x)
        >>> print(k)
        [ 0.          0.25132741  0.50265482 -0.50265482 -0.25132741]
        """
        n = len(x)

        high = (n - 1) // 2
        low = - (n // 2)

        l = n * (x[1] - x[0])

        k = np.concatenate((np.arange(0, high + 1), np.arange(low, 0))) * 2 * np.pi / l

        return k

    def calc_k2(self):
        """Calculates the squared wavenumber.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        numpy.ndarray
            The squared wavenumber.
        """
        return sum([self.k[i] ** 2 for i in range(len(self.k))])

    def calc_Gaussian_filter_f(self, a0=None):
        """Calculate Gaussian filter in Fourier space.

        This method computes a Gaussian filter in Fourier space using the formula:
        exp(-1/2 * a0^2 * k^2) where k^2 is the squared wavenumber.

        Parameters
        ----------
        a0 : float, optional
            Filter width parameter. If None, uses the instance's a0 attribute.

        Returns
        -------
        numpy.ndarray
            The Gaussian filter in Fourier space.

        Notes
        -----
        The filter is calculated using the squared wavenumber obtained from 
        the calc_k2() method.
        """

        if a0 is None:
            a0 = self.a0

        return np.exp(-1 / 2 * a0 ** 2 * self.calc_k2())

    def calc_determinant_field(self, psi: list[np.ndarray, np.ndarray]) -> np.ndarray:
        """Calculate the determinant transformation of a given field.

        The D-field used in the calculation of the topologicel defect density fields. 
        See the ComFiT documentation for more information:
        https://comfitlib.com/TopologicalDefects/
        
        Parameters
        ----------
        psi : list of numpy.ndarray
            A list of two psi fields.
        
        Returns
        -------
        numpy.ndarray
            The determinant field.
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

    def calc_defect_density(self, psi: list[np.ndarray, np.ndarray], psi0=1) -> np.ndarray:
        """Calculate the defect density of a given psi field.
        
        Parameters
        ----------
        psi : list of numpy.ndarray
            A list of two psi fields.
        psi0 : float, optional
            The value of psi_0. Default is 1.
            
        Returns
        -------
        numpy.ndarray
            The defect density of the psi field.
        """

        return 1 / (np.pi * psi0 ** 2) * self.calc_determinant_field(psi)

    def calc_defect_density_singular(self, psi: np.ndarray, psi0=1) -> np.ndarray:
        """Calculate the singular defect density for a given psi field.

        Parameters
        ----------
        psi : numpy.ndarray
            The field psi.
        psi0 : float, optional
            The reference value of psi. Default is 1.

        Returns
        -------
        numpy.ndarray
            The defect density for the given psi value.
        """
        return self.calc_defect_density(psi, 1) * self.calc_delta_function(psi, psi0)

    def calc_defect_velocity_field(
        self,
        psi: list[np.ndarray, np.ndarray],
        dt_psi: list[np.ndarray, np.ndarray]
    ) -> list[np.ndarray, np.ndarray]:
        """Calculate the velocity field of the defects in the psi field.

        Parameters
        ----------
        psi : list of numpy.ndarray
            The psi field.
        dt_psi : list of numpy.ndarray
            The time derivative of the psi field.

        Returns
        -------
        list of numpy.ndarray
            The velocity field of the defects.
        """
        if self.dim == 2:
            if len(psi) == 2:
                # Input to exclude region
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

                # print(region_to_ignore.shape)
                # print(Vx.shape)
                # print(Vy.shape)

                Vx[region_to_ignore] = 0
                Vy[region_to_ignore] = 0

                return [Vx, Vy]

        elif self.dim == 3:
            if len(psi) == 2:
                # Input to exclude region
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
        """Calculates the conserved current of the superfluid density.
        
        Parameters
        ----------
        psi : numpy.ndarray
            The vector field that we find the density of.
        dt_psi : numpy.ndarray
            The time derivative of psi.
        psi_0 : float or numpy.ndarray, optional
            The equilibrium state. Default is 0.
        
        Returns
        -------
        list of numpy.ndarray
            Components of the conserved current.
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

    def calc_delta_function(self, psi: list[np.ndarray, np.ndarray], psi0=1):
        """Calculate the delta function for a given wavefunction.

        Parameters
        ----------
        psi : list of numpy.ndarray
            The wavefunction.
        psi0 : float, optional
            The width of the wavefunction. Default is 1.

        Returns
        -------
        numpy.ndarray
            The value of the delta function.
        """
        width = psi0 / 2
        n = len(psi)
        if self.dim == 2:
            if n == 2:
                psi2 = psi[0] ** 2 + psi[1] ** 2
                return 1 / (2 * np.pi * width ** 2) * np.exp(-psi2 / (2 * width ** 2))

    def calc_region_interval(self, a: float, b: float) -> np.ndarray:
        """Calculate a boolean array indicating whether a point is within an interval.
        
        Parameters
        ----------
        a : float
            The lower bound of the interval.
        b : float
            The upper bound of the interval.
        
        Returns
        -------
        numpy.ndarray
            A boolean array indicating whether a point is within the interval.
            
        Raises
        ------
        Exception
            If the lower bound is greater than the upper bound or if the dimension of the system is not 1.
        """
        if not (a <= b):
            raise Exception("The lower bound must be less than or equal to the upper bound.")
        
        if not (self.dim == 1):
            raise Exception("This function is only valid for 1D systems.")
        
        return (a <= self.x) & (self.x <= b)
        

    def calc_region_disk(self, position: list[float], radius: float) -> np.ndarray:
        """Calculates a boolean array indicating whether a point is within a disk of a given radius.
        
        Parameters
        ----------
        position : list[float]
            The position of the disk
        radius : float
            The radius of the disk
        
        Returns
        -------
        numpy.ndarray
            A boolean array indicating whether a point is within the disk.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 2.
        """
        if self.dim == 2:
            return self.calc_distance_squared_to_point(position) <= radius ** 2

        else:
            raise Exception("Not valid for other dimensions.")

    def calc_region_ball(self, 
            position: list[float], 
            radius: float) -> np.ndarray:
        """Calculate a boolean array indicating whether a point is within a ball of a given radius.
        
        Parameters
        ----------
        position : list[float]
            The position of the ball
        radius : float
            The radius of the ball
        
        Returns
        -------
        numpy.ndarray
            A boolean array indicating whether a point is within the ball.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 3.
        """
        if self.dim == 3:
            return self.calc_distance_squared_to_point(position) <= radius ** 2

        else:
            raise Exception("Not valid for other dimensions.")

    def calc_Gaussian(self, 
            position=None, 
            width=None, 
            top=None,
            value=None):
        """Calculate the Gaussian function.
        
        Parameters
        ----------
        position : list or array, optional
            The center position of the Gaussian. If None, uses the middle of the grid.
        width : float, optional
            The width of the Gaussian. If None, uses the system's a0 parameter.
        top : float, optional
            The maximum value of the Gaussian function.
        value : float, optional
            The integrated value of the Gaussian function.
            If neither top nor value is provided, the function will be normalized to 1.
        
        Returns
        -------
        numpy.ndarray
            The Gaussian function evaluated on the grid.
        
        Notes
        -----
        If top is provided, it scales the height of the Gaussian.
        If value is provided, it scales the integral of the Gaussian.
        """

        if position is None:
            position = self.rmid
        if width is None:
            width = self.a0

        r2 = self.calc_distance_squared_to_point(position)

        if top is not None:
            return top * np.exp(-r2 / (2 * width ** 2))

        if value is None:
            value =1

        return value*(2*np.pi*width**2)**(-self.dim/2)*np.exp(-r2/(2*width**2))

    def calc_distance_squared_to_point(self, position: Union[float, list[float]]) -> np.ndarray:
        """Calculate the squared distance from each point in the grid to a specified position.

        Parameters
        ----------
        position : Union[float, list[float]]
            The target position. Should be a single float for 1D systems or a list of floats for higher dimensions.

        Returns
        -------
        np.ndarray
            An array of the same dimensionality as the system containing the squared distance 
            from each grid point to the specified position. The shape matches the grid resolution.

        Notes
        -----
        This method accounts for periodic boundary conditions by calculating distances in each
        dimension considering the original position, and positions shifted by +/- the system size.
        It then takes the minimum of these distances to find the shortest path under periodicity.
        """

        if self.dim == 1:
            position = [position]
        
        rx2m = (self.x - position[0] - self.size_x) ** 2
        rx2 = (self.x - position[0]) ** 2
        rx2p = (self.x - position[0] + self.size_x) ** 2

        r2 = np.min(np.stack((rx2m, rx2, rx2p)), axis=0).reshape((self.xRes))

        if self.dim > 1:
            r2 = r2.reshape((self.xRes, 1))

            ry2m = (self.y - position[1] - self.size_y) ** 2
            ry2 = (self.y - position[1]) ** 2
            ry2p = (self.y - position[1] + self.size_y) ** 2
            ry2 = np.min(np.stack((ry2m, ry2, ry2p)), axis=0).reshape((1, self.yRes))

            r2 = r2 + ry2
        
        if self.dim > 2:
            r2 = r2.reshape((self.xRes, self.yRes, 1))
            
            rz2m = (self.z - position[2] - self.size_z) ** 2
            rz2 = (self.z - position[2]) ** 2
            rz2p = (self.z - position[2] + self.size_z) ** 2
            rz2 = np.min(np.stack((rz2m, rz2, rz2p)), axis=0).reshape((1, 1, self.zRes))

            r2 = r2 + rz2

        return r2

    def calc_region_cylinder(
        self,
        position: list[float],
        radius: float,
        normal_vector: np.ndarray,
        height: float
    ) -> np.ndarray:
        """Calculates a boolean array indicating whether a point is within a cylinder of a given radius and height.
        
        Parameters
        ----------
        position : list[float]
            The position of the cylinder
        radius : float
            The radius of the cylinder
        normal_vector : np.ndarray
            The normal vector of the cylinder
        height : float
            The height of the cylinder
        
        Returns
        -------
        np.ndarray
            A boolean array indicating whether a point is within the cylinder.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 3.
        """

        if self.dim == 3:
            t = normal_vector / np.linalg.norm(np.array(normal_vector))

            rx = (self.x - position[0])
            # Points that are further away than half the size of the box are moved to the other side
            rx[rx > self.size_x / 2] = rx[rx > self.size_x / 2] - self.size_x
            rx[rx < -self.size_x / 2] = rx[rx < -self.size_x / 2] + self.size_x

            ry = (self.y - position[1])
            ry[ry > self.size_y / 2] = ry[ry > self.size_y / 2] - self.size_y
            ry[ry < -self.size_y / 2] = ry[ry < -self.size_y / 2] + self.size_y

            rz = (self.z - position[2])
            rz[rz > self.size_z / 2] = rz[rz > self.size_z / 2] - self.size_z
            rz[rz < -self.size_z / 2] = rz[rz < -self.size_z / 2] + self.size_z

            zt = rx * t[0] + ry * t[1] + rz * t[2]

            # Project to perpendicular plane vector
            Rt2 = (rx - zt * t[0]) ** 2 + (ry - zt * t[1]) ** 2 + (rz - zt * t[2]) ** 2

            return (zt ** 2 <= height ** 2) & (Rt2 <= radius ** 2)

        else:
            raise Exception("Not valid for other dimensions.")

    def calc_integrate_field(self, field: np.ndarray, region: Optional[int] = None) -> float:
        """Calculates the integrated field value within a specified region.
        
        Parameters
        ----------
        field : numpy.ndarray
            The field array to integrate.
        region : Optional[int], default=None
            If provided, specifies the region to integrate over.
            If None, integrates over the entire field.
        
        Returns
        -------
        float
            The integrated field value within the region.
            
        Raises
        ------
        Exception
            If the dimension of the field is not compatible.
        """

        if region is None:
            return np.sum(field) * self.dV
        else:
            return np.sum(field[region]) * self.dV

    def calc_integrating_factors_f_and_solver(self, omega_f, method: Literal["ETD2RK", "ETD4RK"]) -> tuple:
        """Calculates the integrating factors and the solver for the evolution equation.
        
        Parameters
        ----------
        omega_f : numpy.ndarray
            The value of omega_f in Fourier space.
        method : {'ETD2RK', 'ETD4RK'}
            The method used for evolution.
        
        Returns
        -------
        tuple
            A tuple containing: integrating_factors_f : list (The integrating 
            factors for the selected method), solver : function (The 
            corresponding solver function for the evolution equation).
                
        Raises
        ------
        Exception
            If the specified method is not implemented.
        """
        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop

        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise Exception("This method is not implemented.")

        return integrating_factors_f, solver

    def calc_advect_field(
        self,
        field: np.ndarray,
        u: np.ndarray,
        field_f: Optional[np.ndarray] = None,
        order: int = 3
    ) -> np.ndarray:
        """Advects field according to the provided displacement field u.
        
        Uses a Taylor expansion up to specified order to compute the advection.
        
        Parameters
        ----------
        field : np.ndarray
            The field to be advected.
        u : np.ndarray
            Displacement field.
        field_f : np.ndarray, optional
            Fourier transform of the field. If None, it will be calculated.
        order : int, default=3
            Order of the Taylor expansion. Must be less than or equal to 3.
            
        Returns
        -------
        np.ndarray
            The field after advection by u.
            
        Raises
        ------
        ValueError
            If the order of the Taylor expansion is greater than 3.
        """

        if order > 3:
            raise ValueError("The order of the taylor expansion must be less than or equal to 3.")

        if field_f is None:
            field_f = sp.fft.fftn(field)

        if order > 0:
            for i in range(self.dim):
                # Calculate the derivative
                difield = sp.fft.ifftn(self.dif[i]*field_f)
                # Advect the PFC
                field = field - u[i]*difield

        if order > 1:
            for i in range(self.dim):
                for j in range(i, self.dim):
                    # Calculate the derivative
                    dijfield = sp.fft.ifftn(self.dif[i]*self.dif[j]*field_f)
                    # Advect the PFC
                    field = field + tool_multinom(i,j)*u[i]*u[j]*dijfield
        
        if order > 2:
            for i in range(self.dim):
                for j in range(i, self.dim):
                    for k in range(j, self.dim):
                        # Calculate the derivative
                        dijkfield = sp.fft.ifftn(self.dif[i]*self.dif[j]*self.dif[k]*field_f)
                        # Advect the PFC
                        field = field - tool_multinom(i,j,k)*u[i]*u[j]*u[k]*dijkfield

        return field

    # def calc_advect_field_alot(
    #     self,
    #     field: np.ndarray,
    #     u: np.ndarray,
    #     field_f: Optional[np.ndarray] = None,
    #     order: int = 3
    # ) -> np.ndarray:

    # if self.dim == 2:
        
    #     Ux = u[0]/self.dx
    #     Uy = u[1]/self.dy



    def calc_evolution_integrating_factors_ETD2RK(self, omega_f: np.ndarray, tol: float = 10 ** (-5)) -> list:
        """Calculates integrating factors for ETD2RK.
        
        Parameters
        ----------
        omega_f : numpy.ndarray
            The value of omega_f in Fourier space.
        tol : float, optional
            Tolerance for when to expand the integrating factors that divide by omega.
            Default is 10^-5.
        
        Returns
        -------
        list
            The list of integrating factors for the ETD2RK method.
        """
        integrating_factors_f = [0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = np.divide(If1 - 1, omega_f, 
                                    out=np.full_like(omega_f, self.dt), #Small omega_f * self.dt values are set to low value limit
                                    where=np.abs(omega_f) * self.dt >= tol)

        integrating_factors_f[2] = np.divide(If1 - 1 - omega_f * self.dt, self.dt * omega_f ** 2,
                                    out=np.full_like(omega_f, self.dt / 2), #Small omega_f * self.dt values are set to low value limit
                                    where=np.abs(omega_f) * self.dt >= tol)

        return integrating_factors_f

    def calc_evolution_integrating_factors_ETD4RK(self, omega_f: np.ndarray, tol: float = 10 ** (-5)) -> list:
        """Calculate the evolution integrating factors using the ETDRK4 method.
            
        Parameters
        ----------
        omega_f : numpy.ndarray
            The value of omega_f in Fourier space.
        tol : float, optional
            Tolerance for when to expand the integrating factors that divide by omega.
            Default is 10^-5.
             
        Returns
        -------
        list
            The list of integrating factors for the ETD4RK method.
        """
        integrating_factors_f = [0, 0, 0, 0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt / 2)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = np.divide(If1 - 1, omega_f,
                        out=np.full_like(omega_f, self.dt / 2), #Small omega_f * self.dt values are set to low value limit
                        where=np.abs(omega_f) * self.dt >= tol)

        integrating_factors_f[2] = np.exp(omega_f * self.dt)
        If2 = integrating_factors_f[2]

        integrating_factors_f[3] = np.divide(
            -4 - omega_f * self.dt + If2 * (4 - 3 * omega_f * self.dt + omega_f ** 2 * self.dt ** 2),
            omega_f ** 3 * self.dt ** 2,
            out=np.full_like(omega_f, self.dt / 6), #Small omega_f * self.dt values are set to low value limit
            where=np.abs(omega_f) * self.dt >= tol)

        integrating_factors_f[4] = np.divide(
            2 * (2 + omega_f * self.dt + If2 * (-2 + omega_f * self.dt)),
            omega_f ** 3 * self.dt ** 2,
            out=np.full_like(omega_f, self.dt / 3), #Small omega_f * self.dt values are set to low value limit
            where=np.abs(omega_f) * self.dt >= tol)

        integrating_factors_f[5] = np.divide(
            -4 - 3 * omega_f * self.dt - omega_f ** 2 * self.dt ** 2 + If2 * (4 - omega_f * self.dt),
            omega_f ** 3 * self.dt ** 2,
            out=np.full_like(omega_f, self.dt / 6), #Small omega_f * self.dt values are set to low value limit
            where=np.abs(omega_f) * self.dt >= tol)
        
        return integrating_factors_f

    def calc_defect_nodes(
        self,
        defect_density: np.ndarray,
        charge_tolerance: Optional[float] = None,
        integration_radius: Optional[float] = None
    ):
        """Calculate the positions and charges of defect nodes based on the defect density.
        
        Parameters
        ----------
        defect_density : numpy.ndarray
            The defect density field. A positive scalar field to be integrated.
        charge_tolerance : float, optional
            The minimum charge value to consider as a defect. If None, uses default values based on dimension.
        integration_radius : float, optional
            The radius to use for integrating around potential defects. If None, uses default values based on dimension.
        
        Returns
        -------
        list
            A list of dictionaries representing the defect nodes. Each dictionary contains:
            - 'position_index': The position index of the defect node in the defect density array.
            - 'position': The physical position of the defect node.
            
        Raises
        ------
        Exception
            If the defect density is not a positive scalar field.
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
            def calc_region(position_index,radius):
                return self.calc_region_disk(position=[
                    self.x.flatten()[position_index[0]],
                    self.y.flatten()[position_index[1]]], 
                    radius=radius)

            def calc_position_from_region(defect_density,region_to_integrate,position_index):
                # Calculate the field to integrate
                field_to_integrate = region_to_integrate * defect_density
                # Normalize the field to integrate
                field_to_integrate = field_to_integrate / np.sum(field_to_integrate)

                # Roll the field so that the defect node is at the center
                Rx = round((self.x.flatten()[position_index[0]] - self.xmid) / self.dx)
                field_to_integrate = np.roll(field_to_integrate, -Rx, axis=0)
                Ry = round((self.y.flatten()[position_index[1]] - self.ymid) / self.dy)
                field_to_integrate = np.roll(field_to_integrate, -Ry, axis=1)
                
                x = np.sum(field_to_integrate * self.x) 
                y = np.sum(field_to_integrate * self.y) 

                x = x + Rx*self.dx
                y = y + Ry*self.dy

                return [x,y]

        elif self.dim == 3:
            if charge_tolerance is None:
                charge_tolerance = 0.5*self.a0**2
                # print("charge tolerance",charge_tolerance)
            if integration_radius is None:
                integration_radius = 2*self.a0

            #Auxiliary functions
            def calc_region(position_index,radius):
                return self.calc_region_ball(position=[
                        self.x.flatten()[position_index[0]],
                        self.y.flatten()[position_index[1]],
                        self.z.flatten()[position_index[2]]], 
                        radius=radius)
            
            def calc_position_from_region(defect_density,region_to_integrate,position_index):
                # Calculate the field to integrate
                field_to_integrate = region_to_integrate * defect_density
                # Normalize the field to integrate
                field_to_integrate = field_to_integrate / np.sum(field_to_integrate)

                # Roll the field so that the defect node is at the center
                Rx = round((self.x.flatten()[position_index[0]] - self.xmid) / self.dx)
                field_to_integrate = np.roll(field_to_integrate, -Rx, axis=0)
                Ry = round((self.y.flatten()[position_index[1]] - self.ymid) / self.dy)
                field_to_integrate = np.roll(field_to_integrate, -Ry, axis=1)
                Rz = round((self.z.flatten()[position_index[2]] - self.zmid) / self.dz)
                field_to_integrate = np.roll(field_to_integrate, -Rz, axis=2)

                x = np.sum(field_to_integrate * self.x)
                y = np.sum(field_to_integrate * self.y)
                z = np.sum(field_to_integrate * self.z)

                x = x + Rx*self.dx
                y = y + Ry*self.dy
                z = z + Rz*self.dz

                return [x,y,z]

        #Region to search for defect nodes
        region_to_search = np.ones(self.dims)

        # Calculate the point where defect density is largest
        position_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

        # Integrate the defect density around this point (i.e. in a disk/ball around)
        region_to_integrate = calc_region(position_index,
                                        radius=integration_radius)

        charge = self.calc_integrate_field(defect_density, region_to_integrate)

        # print("Charge: ", charge)
        # print("Charge tolerance: ", charge_tolerance)

        while charge > charge_tolerance:
            defect_node = {}
            defect_node['position_index'] = position_index

            defect_node['position'] = calc_position_from_region(defect_density,region_to_integrate,position_index)

            # print("Defect node position: ", defect_node['position'])
            defect_nodes.append(defect_node)

            region_to_exclude_from_search = calc_region(position_index,
                                        radius=2*integration_radius)
        
            region_to_search[region_to_exclude_from_search] = 0

            position_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

            region_to_integrate = calc_region(position_index,
                                        radius=integration_radius)
            
            charge = self.calc_integrate_field(defect_density, region_to_integrate)

        return defect_nodes