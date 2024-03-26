import numpy as np
import scipy as sp

from comfit.tools.tool_math_functions import tool_multinom

class BaseSystemCalc:
    # CALCULATION FUNCTIONS

    # Calculation of angle fields for vortices of different types
    def calc_angle_field_single_vortex(self,
                                       position=None,
                                       charge=1):
        """
        Calculate the angle field due to a single vortex.

        Input:
            position (list, optional): The position of the vortex. Defaults to None.
            charge (int, optional): The charge of the vortex. Defaults to 1.
        
        Output:
            numpy.ndarray: The angle field calculated for the vortex.
            
        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]
        
        theta = charge * np.arctan2(self.y - position[1], self.x - position[0])

        return np.mod(theta + np.pi, 2 * np.pi) - np.pi

    def calc_angle_field_vortex_dipole(self,
                                       dipole_vector=None,
                                       dipole_position=None):
        """
        Calculates the angle field for a double vortex system.
        Input:
            dipole_vector (list, optional): The dipole vector. Defaults to None.
            dipole_position (list, optional): The position the center of mass of the dipole. Defaults to None.
        Raises:
            Exception: If the dimension of the system is not 2.
        Output:
            np.ndarray: The calculated angle field for the double vortex system.
        """

        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if dipole_vector is None:
            dipole_vector = [self.xmax / 3, 0]

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
        width = 0.2 * np.min([self.xmax, self.ymax])
        radius = 0.4 * np.min([self.xmax, self.ymax])

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

    def calc_angle_field_vortex_ring(self, position=None, radius=None, normal_vector=[0, 0, 1]):
        """
        Calculates the angle field for a vortex ring.
        Input:
            position (list, optional): The position of the vortex ring. Defaults to None.
            radius (float, optional): The radius of the vortex ring. Defaults to None.
            normal_vector (list, optional): The normal vector of the vortex ring. Defaults to [0,0,1].
        Output:
            numpy.ndarray: The calculated angle field for the vortex ring.
        """
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

        Input:
        - x : 1D array of x-positions.

        Output:
            numpy array: 1D array of wavenumbers with all the modes for the given x-array,
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

    def calc_Gaussian_filter_f(self, a0=None):

        if a0 is None:
            a0 = self.a0

        return np.exp(-1 / 2 * a0 ** 2 * self.calc_k2())

    def calc_determinant_field(self, psi):
        """
        Calculate the determinant transformation of a given field

        Input:
            psi (list): A list of two psi fields.

        Output:
            numpy.ndarray: The defect density of the psi field.
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

        Input:
            psi (list): A list of two psi fields.
            psi0 (float, optional): The value of psi_0. Defaults to 1.

        Output:
            np.ndarray: The defect density of the psi field.
        """

        return 1 / (np.pi * psi0 ** 2) * self.calc_determinant_field(psi)

    def calc_defect_density_singular(self, psi, psi0=1):
        """
        Calculate the singular defect density for a given psi field.

        Input:
            psi (float): The value of psi.
            psi0 (float, optional): The reference value of psi. Defaults to 1.
        Output:
            np.ndarray: The defect density for the given psi value.
        """
        return self.calc_defect_density(psi, 1) * self.calc_delta_function(psi, psi0)

    def calc_defect_velocity_field(self, psi, dt_psi):
        """
        Calculates the velocity field of the defects in the psi field.

        Input:
            psi: The psi field
            dt_psi: The time derivative of the psi field

        Output:
            np.ndarray: The velocity field of the defects
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

                print(region_to_ignore.shape)
                print(Vx.shape)
                print(Vy.shape)

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
        """
        Calculates the conserved current of the superfluid density

        Input:
            psi (numpy.ndarray) the vector field that we find the density of
            dt_psi (numpy.ndarray) the time derivative of psi
            psi_0 (floar or numpy.ndarray, optional) the equilibrium state
        
        Output:
            np.ndarray: Components of the conserved current
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

        Input:
            psi (list): The wavefunction.
            psi0 (float): The width of the wavefunction. Default is 1.

        Output:
            np.ndarray: The value of the delta function.
        """
        width = psi0 / 2
        n = len(psi)
        if self.dim == 2:
            if n == 2:
                psi2 = psi[0] ** 2 + psi[1] ** 2
                return 1 / (2 * np.pi * width ** 2) * np.exp(-psi2 / (2 * width ** 2))

    def calc_region_interval(self,a,b):
        """
        Calculates a boolean array indicating whether a point is within an interval.
        
        Input:
            a: The lower bound of the interval
            b: The upper bound of the interval
        
        Output:
            np.ndarray: A boolean array indicating whether a point is within the interval
        """
        if not (a <= b):
            raise Exception("The lower bound must be less than or equal to the upper bound.")
        
        if not (self.dim == 1):
            raise Exception("This function is only valid for 1D systems.")
        
        return (a <= self.x) & (self.x <= b)
        

    def calc_region_disk(self, position, radius):
        """
        Calculates a boolean array indicating whether a point is within a disk of a given radius.
        
        Input:
            position: The position of the disk
            radius: The radius of the disk
        
        Output:
            np.ndarray: A boolean array indicating whether a point is within the disk
        """
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
        
        Input:
            position: The position of the ball
            radius: The radius of the ball
        
        Output:
            np.ndarray: A boolean array indicating whether a point is within the ball
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

    def calc_Gaussian(self, **kwargs):
        """
        Calculated the Gaussian function 
        Input:
            kwargs:
            - position
            - width
            - top (the top of the Gaussian function)
            - value (the value of the integrated Gaussian function)
            If neither top nor value is provided, the function will be normalized to 1
        """
        position = kwargs.get('position',self.rmid)
        width = kwargs.get('width',self.a0)

        r2 = self.calc_distance_squared_to_point(position)
        
        if 'top' in kwargs:
            return kwargs['top']*np.exp(-r2/(2*width**2))
        else:
            value = kwargs.get('value',1)
            return value*(2*np.pi*width**2)**(-self.dim/2)*np.exp(-r2/(2*width**2))

    def calc_distance_squared_to_point(self,position):
        """
        Calculates the distance to a point given
        """

        if self.dim == 1:
            position = [position]
        
        delta_x = self.xmax-self.xmin
        rx2m = (self.x - position[0] - delta_x) ** 2
        rx2 = (self.x - position[0]) ** 2
        rx2p = (self.x - position[0] + delta_x) ** 2

        r2 = np.min(np.stack((rx2m, rx2, rx2p)), axis=0).reshape((self.xRes))

        if self.dim > 1:
            r2 = r2.reshape((self.xRes, 1))

            delta_y = self.ymax-self.ymin
            ry2m = (self.y - position[1] - delta_y) ** 2
            ry2 = (self.y - position[1]) ** 2
            ry2p = (self.y - position[1] + delta_y) ** 2
            ry2 = np.min(np.stack((ry2m, ry2, ry2p)), axis=0).reshape((1, self.yRes))

            r2 = r2 + ry2
        
        if self.dim > 2:
            r2 = r2.reshape((self.xRes, self.yRes, 1))
            
            delta_z = self.zmax-self.zmin
            rz2m = (self.z - position[2] - delta_z) ** 2
            rz2 = (self.z - position[2]) ** 2
            rz2p = (self.z - position[2] + delta_z) ** 2
            rz2 = np.min(np.stack((rz2m, rz2, rz2p)), axis=0).reshape((1, 1, self.zRes))

            r2 = r2 + rz2

        return r2

    def calc_region_cylinder(self, position, radius, normal_vector, height):
        """
        Calculates a boolean array indicating whether a point is within a cylinder of a given radius and height.
        
        Input:
            position: The position of the cylinder
            radius: The radius of the cylinder
            normal_vector: The normal vector of the cylinder
            height: The height of the cylinder
        
        Output:
            np.ndarray: A boolean array indicating whether a point is within the cylinder
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

        Input:
            field (numpy.ndarray): The field array.
            index (tuple, optional): The indices of the center point in the field. Defaults to None.
            radius (float, optional): The radius of the region. Defaults to None.

        Input:
            tuple or float: If index is provided, returns a tuple containing the integrated field value
                            within the region and a boolean array indicating the region. If index is None,
                            returns the integrated field value within the entire field.

        Output:
            float: The integrated field value within the region.
                                            
        Raises:
            Exception: If the dimension of the field is not 2.
        """

        if region is None:
            return np.sum(field) * self.dV
        else:
            return np.sum(field[region]) * self.dV

    def calc_integrating_factors_f_and_solver(self, omega_f, method):
        """
        Calculates the integrating factors and the solver for the evolution equation.
        
        Input:
            omega_f: The value of omega_f
            method: The method used for evolution
        
        Output:
            integrating_factors_f: The integrating factors
            solver: The solver for the evolution equation
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

    def calc_advect_field(self, field, u, field_f = None, order = 3):
        """
        Advects field accodring to the provided displacement field u, with a taylor expansion up to order 3.

        Input: 
            u: Displacement field (numpy.ndarray)
            taylor_order: Order of the taylor expansion

        Output:
            np.ndarray: The field after advection by u
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


    def calc_evolution_integrating_factors_ETD2RK(self, omega_f, tol=10 ** (-4)):
        """
        Calculates integrating factors for ETD2RK
        
        Input:
            omega_f (numpy.ndarray): the value of omega_f
            tol (float, optional): tolerance for when to expand the integrating factors that divide by omega
        
        Output:
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

        Input:
            omega_f (numpy.ndarray): The value of omega_f.
            tol (float,optional): tolerance for when to expand the integrating factors that divide by omega
        
        Output:
            list: The list of integrating factors.
        """
        integrating_factors_f = [0, 0, 0, 0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt / 2)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = (If1 - 1) / omega_f

        integrating_factors_f[2] = np.exp(omega_f * self.dt)
        If2 = integrating_factors_f[2]

        integrating_factors_f[3] = 1 / ( omega_f ** 3 * self.dt ** 2) \
                                   * (-4 - omega_f * self.dt + If2 * (
                4 - 3 * omega_f * self.dt +  omega_f ** 2 * self.dt ** 2))

        integrating_factors_f[4] = 2 / (omega_f ** 3 * self.dt ** 2) \
                                   * (2 + omega_f*self.dt + If2 * (-2 + omega_f * self.dt))

        integrating_factors_f[5] = 1 / (omega_f ** 3 * self.dt ** 2) \
                                   * (-4 - 3 * omega_f * self.dt - omega_f ** 2 * self.dt ** 2 + If2 * (
                4 - omega_f*self.dt))


        #Small omega_f limits
        integrating_factors_f[1][np.abs(omega_f) < tol] = self.dt / 2
        integrating_factors_f[3][np.abs(omega_f) < tol] = self.dt / 6
        integrating_factors_f[4][np.abs(omega_f) < tol] = self.dt / 3
        integrating_factors_f[5][np.abs(omega_f) < tol] = self.dt / 6

        return integrating_factors_f

    def calc_defect_nodes(self, defect_density,
                          charge_tolerance=None,
                          integration_radius=None):
        """
        Calculate the positions and charges of defect nodes based on the defect density.
        
        Input:
            defect_density (numpy.ndarray): The defect density field. A positive scalar field to be integrated.
        
        Output:
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
        position_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

        # Integrate the defect density around this point (i.e. in a disk/ball around)
        region_to_integrate = calc_region(position_index,
                                        radius=integration_radius)

        charge = self.calc_integrate_field(defect_density, region_to_integrate)

        while charge > charge_tolerance:
            defect_node = {}
            defect_node['position_index'] = position_index

            defect_node['position'] = calc_position_from_region(defect_density,region_to_integrate)

            defect_nodes.append(defect_node)

            region_to_exclude_from_search = calc_region(position_index,
                                        radius=2*integration_radius)
        
            region_to_search[region_to_exclude_from_search] = 0

            position_index = np.unravel_index(np.argmax(defect_density*region_to_search), defect_density.shape)

            region_to_integrate = calc_region(position_index,
                                        radius=integration_radius)
            
            charge = self.calc_integrate_field(defect_density, region_to_integrate)

        return defect_nodes