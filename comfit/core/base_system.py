import numpy as np
import matplotlib.pyplot as plt
from comfit.tools.tool_colormap_angle import tool_colormap_angle


class BaseSystem:
    def __init__(self, dimension, xRes=101, dx=1.0, yRes=101, dy=1.0, zRes=101, dz=1.0, dt=0.1):
        self.dim = dimension
        self.xRes = xRes
        self.yRes = 1
        self.zRes = 1

        if dimension > 1:
            self.yRes = yRes

        if dimension > 2:
            self.zRes = zRes

        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.dt = dt

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
        self.xmidi = (1 + self.xRes) // 2 - 1
        self.xmid = self.x[self.xmidi]

        self.ymidi = (1 + self.yRes) // 2 - 1
        self.ymid = self.y[self.ymidi]

        self.zmidi = (1 + self.zRes) // 2 - 1
        self.zmid = self.z[self.zmidi]

        self.midi = self.xRes * self.yRes * (self.zmidi - 1) + self.yRes * (self.xmidi - 1) + self.ymidi
        self.rmid = [self.xmid, self.ymid, self.zmid]

        # Max positions
        self.xmax = self.x[-1] + self.dx
        self.ymax = self.y[-1] + self.dy if self.dim > 1 else 0
        self.zmax = self.z[-1] + self.dz if self.dim > 2 else 0

        # Fourier modes
        self.k = [self.calc_wavenums(self.x), 0, 0]
        if self.dim > 1:
            self.k[1] = self.calc_wavenums(self.y)
        if self.dim > 2:
            self.k[2] = self.calc_wavenums(self.z)

        # Derivatives
        self.dif = [1j * ki for ki in self.k]

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

    # Calculation of angle fields for vortices of different types
    def calc_angle_field_single_vortex(self,
                                       position=None,
                                       charge=1):
        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]

        x = self.y.reshape((self.xRes, 1))
        y = self.y.reshape((1, self.yRes))

        theta = charge*np.arctan2(y - position[1],x - position[0])

        return theta

    def calc_angle_field_double_vortex(self,
                                       position1=None,
                                       position2=None):

        if self.dim != 2:
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position1 is None:
            position1 = [self.xmax/3, self.ymid]

        if position2 is None:
            position2 = [2*self.xmax/3,self.ymid]

        theta1 = self.calc_angle_field_single_vortex(position1)
        theta2 = self.calc_angle_field_single_vortex(position2,charge=-1)

        return np.mod(theta1+theta2+np.pi,2*np.pi)-np.pi

    def calc_wavenums(self,x):
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

        high = n // 2
        low = - (n - 1) // 2

        l = n * (x[1] - x[0])

        k = np.concatenate((np.arange(0, high + 1), np.arange(low, 0))) * 2 * np.pi / l

        if n % 2 == 0:
            k[n // 2] = -k[n // 2]

        return k

    def calc_k2(self):
        if self.dim == 1:
            k2 = self.k[0] ** 2
        elif self.dim == 2:
            k2 = (self.k[0] ** 2).reshape(self.xRes, 1) + \
                 (self.k[1] ** 2).reshape(1, self.yRes)
        elif self.dim == 3:
            k2 = (self.k[0] ** 2).reshape(self.xRes, 1, 1) + \
                 (self.k[1] ** 2).reshape(1, self.yRes, 1) + \
                 (self.k[2] ** 2).reshape(1, 1, self.zRes)
        return k2




    #plotting functions
    def plot_angle_field(self,field):
        X,Y = np.meshgrid(self.x,self.y,indexing='ij')

        custom_colormap = tool_colormap_angle()

        mesh = plt.pcolormesh(X,Y,field,shading='auto',cmap=custom_colormap)
        cbar = plt.colorbar(mesh)  # To add a colorbar on the side
        cbar.set_ticks(np.array([-np.pi,-2*np.pi/3,-np.pi/3,0,np.pi/3,2*np.pi/3,np.pi]))
        cbar.set_ticklabels([r'$-\pi$',r'$-2\pi/3$',r'$-\pi/3$',r'$0$',r'$\pi/3$',r'$2\pi/3$',r'$\pi$'])
        plt.title("Angle field")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.show()

    def plot_field(self,field,ax):

        if self.dim == 2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            surf = ax.plot_surface(X, Y, field)

    def plot_complex_field(self,complex_field,ax):



        if self.dim==2:
            rho = np.abs(complex_field)
            theta = np.angle(complex_field)

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            custom_colormap = tool_colormap_angle()

            # Get the colors from a colormap (e.g., hsv, but you can choose any other)
            colors = plt.cm.hsv((theta + np.pi) / (2 * np.pi))  # Normalizing theta to [0, 1]

            surf = ax.plot_surface(X, Y, rho, facecolors=colors)

            # mappable = plt.cm.ScalarMappable(cmap=custom_colormap)
            # mappable.set_array([])
            # mappable.set_clim(-np.pi, np.pi)
            # cbar = plt.colorbar(mappable, ax=ax)
            # cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
            # cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

            # plt.title("Angle field")
            # plt.xlabel("X-axis")
            # plt.ylabel("Y-axis")



        else:
            raise Exception("This plotting function not yet configured for other dimension")


    # Time evolution function
    def evolve_ETDRK2_loop(self,integrating_factors_f,field,field_f,number_of_pred_it_steps=2):

        N0_f = self.calc_nonlinear_evolution_term_f(field)
        #This needs to be generalized

        for i in range(number_of_pred_it_steps):
            if i==1:
                dN_f = 0
            else:
                dN_f = self.calc_nonlinear_evolution_term_f(field) - N0_f

            #print(N0_f)
            field_f_pred = integrating_factors_f[0]*field_f +\
                integrating_factors_f[1]*N0_f +\
                integrating_factors_f[2]*dN_f

            #Can this be simplified? (Vidar 08.09.23)
            if self.dim==1:
                field_f_pred[0] = field_f[0]
            elif self.dim==2:
                field_f_pred[0,0] = field_f[0,0]
            elif self.dim==3:
                field_f_pred[0,0,0] = field_f[0,0,0]

            field = np.fft.ifftn(field_f_pred)

        return field, field_f_pred


