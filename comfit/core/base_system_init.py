import numpy as np
import scipy as sp
from comfit.tool.tool_print_in_color import tool_print_in_color
from comfit.tool.tool_configure_axis import tool_configure_axis

class BaseSystemInit:
    """ Initialization methods for the base system class"""
    def __init__(self, dim: int, **kwargs):
        """
        Initialize the class with the given parameters.

        Args:
            dimension (int): The dimension of the system. Must be 1, 2, or 3.
            **kwargs: Additional keyword arguments, see 
            https://comfitlib.com/ClassBaseSystem/

        Returns:
            None
        """

        if dim not in [1, 2, 3]:
            raise ValueError('Dimension must be 1, 2, or 3.')

        self.dim = dim

        # Checking what x-values were provided
        self.xlim = kwargs.get('xlim', None)
        self.xmin = kwargs.get('xmin', None)
        self.xmax = kwargs.get('xmax', None)
        self.xRes = kwargs.get('xRes', None)
        self.dx = kwargs.get('dx', None)

        self.xlim, self.xmin, self.xmax, self.xRes, self.dx = tool_configure_axis(self.dim, 'x',self.xlim, self.xmin, self.xmax, self.xRes, self.dx)

        # Checking what y-values were provided
        self.ylim = kwargs.get('ylim', None)
        self.ymin = kwargs.get('ymin', None)
        self.ymax = kwargs.get('ymax', None)
        self.yRes = kwargs.get('yRes', None)
        self.dy = kwargs.get('dy', None)

        self.ylim, self.ymin, self.ymax, self.yRes, self.dy = tool_configure_axis(self.dim, 'y',self.ylim, self.ymin, self.ymax, self.yRes, self.dy)

        # Checking what z-values were provided
        self.zlim = kwargs.get('zlim', None)
        self.zmin = kwargs.get('zmin', None)
        self.zmax = kwargs.get('zmax', None)
        self.zRes = kwargs.get('zRes', None)
        self.dz = kwargs.get('dz', None)

        self.zlim, self.zmin, self.zmax, self.zRes, self.dz = tool_configure_axis(self.dim, 'z',self.zlim, self.zmin, self.zmax, self.zRes, self.dz)

        # Minimum and maximum sizes of the domain
        self.size_x = self.xmax - self.xmin
        self.size_y = self.ymax - self.ymin if self.dim > 1 else 1
        self.size_z = self.zmax - self.zmin if self.dim > 2 else 1

        self.volume = self.size_x * self.size_y * self.size_z

        if self.dim == 1:
            self.size_min = self.size_x
            self.size_max = self.size_x
        elif self.dim == 2:
            self.size_min = np.min([self.size_x, self.size_y])
            self.size_max = np.max([self.size_x, self.size_y])
        elif self.dim == 3:
            self.size_min = np.min([self.size_x, self.size_y, self.size_z])
            self.size_max = np.max([self.size_x, self.size_y, self.size_z])
        

        #  Setting default value of time
        self.time = kwargs.get('time', 0)
        self.dt = kwargs.get('dt', 0.1)

        # Construct parameters
        self.x = np.linspace(self.xmin, self.xmax-self.dx, self.xRes)
        self.y = np.linspace(self.ymin, self.ymax-self.dy, self.yRes)
        self.z = np.linspace(self.zmin, self.zmax-self.dz, self.zRes)

        self.Res = self.xRes * self.yRes * self.zRes
        if self.dim == 1:
            self.dims = self.xRes
        elif self.dim == 2:
            self.dims = [self.xRes, self.yRes]
        elif self.dim == 3:
            self.dims = [self.xRes, self.yRes, self.zRes]

        self.a0 = 1  # System length scale

        # Helpful midpoints and their indices
        self.xmidi = (self.xRes) // 2 
        self.xmid = self.x[self.xmidi]

        self.ymidi = (self.yRes) // 2 
        self.ymid = self.y[self.ymidi]

        self.zmidi = (self.zRes) // 2 
        self.zmid = self.z[self.zmidi]

        self.midi = self.xRes * self.yRes * (self.zmidi - 1) + self.yRes * (self.xmidi - 1) + self.ymidi
        if self.dim == 1:
            self.rmid = self.xmid
            self.zero_index = 0 
        elif self.dim == 2:
            self.rmid = [self.xmid, self.ymid]
            self.zero_index = (0,0)
        elif self.dim == 3:
            self.rmid = [self.xmid, self.ymid, self.zmid]
            self.zero_index = (0,0,0)

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


        # Plot lib
        self.plot_lib = kwargs.get('plot_lib', 'plotly')

    def __str__(self) -> str:
        """Return a string representation of the class.
        
        Input:
            None

        Returns:
            A string representation of the class.
        """
        description = f"BaseSystem instance\n"
        # Start with the dimension of the system
        description = f"System Dimension: {self.dim}\n"

        # Add x-axis properties
        description += f"X-Axis Limits: [{self.xmin}, {self.xmax}], Resolution: {self.xRes}, Delta: {self.dx}\n"

        # Add y-axis properties if dim > 1
        if self.dim > 1:
            description += f"Y-Axis Limits: [{self.ymin}, {self.ymax}], Resolution: {self.yRes}, Delta: {self.dy}\n"

        # Add z-axis properties if dim > 2
        if self.dim > 2:
            description += f"Z-Axis Limits: [{self.zmin}, {self.zmax}], Resolution: {self.zRes}, Delta: {self.dz}\n"

        # Add time properties
        description += f"Current Time: {self.time}, Time Step: {self.dt}"

        return description

