import numpy as np
import scipy as sp


class BaseSystemInit:
    def __init__(self, dim, **kwargs):
        """
        Initialize the class with the given parameters.

        Input:
            dimension (int): The dimension of the system. Must be 1, 2, or 3.
            **kwargs: Additional keyword arguments, see 
            https://vidarsko.github.io/ComFiT/ClassBaseSystem/

        Output:
            None
        """
        if dim not in [1, 2, 3]:
            raise ValueError('Dimension must be 1, 2, or 3.')

        self.dim = dim

        if 'xlim' in kwargs:
            # Providing xlim trumps providing xmin and xmax
            self.xmin = kwargs['xlim'][0]
            self.xmax = kwargs['xlim'][1]
            
            self.xRes = kwargs.get('xRes', 101)
            self.dx = (self.xmax - self.xmin) / self.xRes
        
        else:
            self.xmin = kwargs.get('xmin', 0)
            self.xRes = kwargs.get('xRes', 101)

            self.dx = kwargs.get('dx', 1.0)

            if 'xmax' in kwargs:
                self.xmax = kwargs['xmax']
            else:
                self.xmax = self.xmin + self.xRes * self.dx

        self.xlim = [self.xmin, self.xmax]

        # Providing dx trumps providing xRes
        if 'dx' in kwargs:
            self.dx = kwargs['dx']
            self.xRes = round((self.xmax - self.xmin) / self.dx)

        # Setting default values for y and z for 1 dimensional systems
        self.ymin = kwargs.get('ymin', 0)
        self.zmin = kwargs.get('zmin', 0)
        
        self.dy = kwargs.get('dy', 1.0)
        self.dz = kwargs.get('dz', 1.0)

        self.ymax = kwargs.get('ymax',1)
        self.zmax = kwargs.get('zmax',1)

        self.yRes = kwargs.get('yRes', 1)
        self.zRes = kwargs.get('zRes', 1)

        if self.dim > 1:
            if 'ylim' in kwargs:
                self.ymin = kwargs['ylim'][0]
                self.ymax = kwargs['ylim'][1]
                
                self.yRes = kwargs.get('yRes', 101)
                self.dy = (self.ymax - self.ymin) / self.yRes
            
            else:
                self.ymin = kwargs.get('ymin', 0)
                self.yRes = kwargs.get('yRes', 101)
                self.dy = kwargs.get('dy', 1.0)

                if 'ymax' in kwargs:
                    self.ymax = kwargs['ymax']
                else:
                    self.ymax = self.ymin + self.yRes * self.dy

            # Providing dy trumps providing yRes
            if 'dy' in kwargs:
                self.dy = kwargs['dy']
                self.yRes = round((self.ymax - self.ymin) / self.dy)
        
        if self.dim > 2:
            if 'zlim' in kwargs:
                self.zmin = kwargs['zlim'][0]
                self.zmax = kwargs['zlim'][1]
                
                self.zRes = kwargs.get('zRes', 101)
                self.dz = (self.zmax - self.zmin) / self.zRes

            else:
                self.zmin = kwargs.get('zmin', 0)
                self.zRes = kwargs.get('zRes', 101)
                self.dz = kwargs.get('dz', 1.0)

                if 'zmax' in kwargs:
                    self.zmax = kwargs['zmax']
                else:
                    self.zmax = self.zmin + self.zRes * self.dz

            # Providing dz trumps providing zRes
            if 'dz' in kwargs:
                self.dz = kwargs['dz']
                self.zRes = round((self.zmax - self.zmin) / self.dz)

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
            self.zero_index = 0 #Nothing comment
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

    def __str__(self):
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

