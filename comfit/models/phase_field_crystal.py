import numpy as np
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from scipy.optimize import fsolve


class PhaseFieldCrystal(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Nothing here yet
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.dislocation_charges = np.array(
            [[np.round(np.dot(an, qn) / (2 * np.pi), decimals=8) for qn in self.q] for an in self.a])

    # CONFIGURATION FUNCTIONS
    def conf_PFC_from_amplitudes(self, eta=None):

        self.psi = self.psi0

        if eta is None:
            eta = self.eta0

        # if self.dim == 1:
        #     X = self.x
        #
        # elif self.dim == 2:
        #     X, Y = np.meshgrid(self.x, self.y, indexing='ij')
        #
        # elif self.dim == 3:
        #     X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        for n in range(self.number_of_reciprocal_modes):

            if self.dim == 1:
                self.psi += 2 * eta[n] * np.exp(1j * self.q[n][0] * self.x)

            elif self.dim == 2:
                self.psi += 2 * eta[n] * np.exp(1j * (self.q[n][0] * self.x + self.q[n][1] * self.y))

            elif self.dim == 3:
                self.psi += 2 * eta[n] * np.exp(1j * (self.q[n][0] * self.x + self.q[n][1] * self.y + self.q[n][2] * self.z))

        self.psi = np.real(self.psi)
        self.psi_f = np.fft.fftn(self.psi)

    # EVOLUTION FUNCTIONS
    def evolve_PFC(self, number_of_steps, method='ETD2RK'):

        omega_f = self.calc_omega_f()

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the PFC'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_f,
                                          self.psi, self.psi_f)

        # These steps seem to be necessary for numerical stability (Vidar 18.12.23)
        self.psi = np.real(self.psi)
        self.psi_f = np.fft.fftn(self.psi)

    def evolve_PFC_hydrodynamic(self, number_of_steps, 
                                method = 'ETD2RK',
                                gamma_S = 2**-2,
                                rho0 = 2**-2):
        """
        We extend psi to also contain the hydrodynamic field in components psi[1],psi[2],psi[3]
        """
        
        if hasattr(self,'velocity_field'):
            print("Velocity field established.")
            pass
        else:
            print("I am initializing a velocity field")
            self.velocity_field = True
            self.psi = np.array([self.psi]+[np.zeros_like(self.psi)]*self.dim)
            self.psi_f = np.array([self.psi_f]+[np.zeros_like(self.psi_f)]*self.dim)
            # print("psi shape", self.psi.shape)

        self.gamma_S = gamma_S
        self.rho0 = rho0

        omega_f = self.calc_omega_hydrodynamic_f()
        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the hydrodynamic PFC'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_hydrodynamic_evolution_function_f,
                                          self.psi, self.psi_f)
            
        self.psi = np.real(self.psi)
        self.psi_f = np.fft.fftn(self.psi, axes = (range ( - self.dim , 0) ) )


    # CALCULATION FUNCTIONS

    def calc_omega_hydrodynamic_f(self):
        return np.array([self.calc_omega_f()]+[-self.gamma_S/self.rho0*np.ones(self.dims)]*self.dim)


    def calc_PFC_free_energy_density_and_chemical_potential(self,field=None,field_f=None):

        if field is None:
            field = self.psi
            field_f = self.psi_f

        psi_f = field_f
        
        # print("field shape",field.shape)
        # print("field_f shape",field_f.shape)

        psi = field 
        psi2 = field**2 
        psi3 = psi2*psi
        psi4 = psi2**2

        free_energy_density = 1/2*self.calc_Lpsi(psi_f)**2 \
            + 1/2*self.r*psi2 + 1/3*self.t*psi3 + 1/4*self.v*psi4
        
        # print("Lpsi shape",self.calc_Lpsi(psi_f).shape)

        chemical_potential = self.calc_L2psi(psi_f)*psi \
            + self.r*psi + self.t*psi2 + self.v*psi3
        
        # print("L2psi shape",self.calc_L2psi(psi_f).shape)
        
        # print("Free energy density shape", free_energy_density.shape)
        # print("Chem pot shape", chemical_potential.shape)
        return free_energy_density, chemical_potential

    def calc_nonlinear_hydrodynamic_evolution_function_f(self, field):

        field_f = np.fft.fftn(field, axes =( range ( - self . dim , 0) ))

        k2 = self.calc_k2()

        N0_f = -k2*np.fft.fftn(self.t * field[0] ** 2 + self.v * field[0] ** 3) \
            - np.fft.fftn(sum([field[i+1]*np.fft.ifftn(self.dif[i]*field_f[0]) for i in range(self.dim)]))
        
        force_density_f = self.calc_stress_divergence_f(field_f[0])

        return np.array([N0_f] + [1/self.rho0*force_density_f[i] for i in range(self.dim)])

    def calc_nonlinear_evolution_function_f(self, psi):
        return -self.calc_k2()*np.fft.fftn(self.t * psi ** 2 + self.v * psi ** 3)

    # Initial configuration methods
    def calc_amplitudes_with_dislocation(self, eta=None, x=None, y=None, dislocation_type=1):

        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a single point dislocation.")

        if x == None:
            x = self.xmid
        if y == None:
            y = self.ymid
        if eta == None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type - 1]
        for n in range(self.number_of_reciprocal_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * self.calc_angle_field_single_vortex([x, y], charge=sn[n]))

        return eta

    def calc_amplitudes_with_dislocation_dipole(self, eta=None,
                                                x1=None, y1=None,
                                                x2=None, y2=None,
                                                dislocation_type=1):
        """
        Inserts a  dislocation dipole in the system corresponding to dislocation type and its negative

        :param eta: The amplitudes currently being manipulated (default: self.eta0)
        :param x1: The x-position of the first dislocation (default: 1/3 xmax)
        :param y1: The y-position of the first dislocation (default: 1/2 ymax)
        :param x2: The x-position of the second dislocation (default: 2/3 xmax)
        :param y2: The y-position of the second dislocation (default: 1/2 ymax)
        :param dislocation_type: the type of dislocation dipole (default: 1)
        :return:
        """
        # TODO: Improve code documentation

        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 to insert a dislocation dipole.")

        if x1 == None:
            x1 = self.xmax / 3
        if y1 == None:
            y1 = self.ymax / 2
        if x2 == None:
            x2 = 2 * self.xmax / 3
        if y2 == None:
            y2 = self.ymax / 2
        if eta == None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type - 1]

        theta = self.calc_angle_field_vortex_dipole(
            dipole_vector=[x2 - x1, y2 - y1],
            dipole_position=[(x1 + x2) / 2, (y1 + y2) / 2])



        for n in range(self.number_of_reciprocal_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * sn[n] * theta)

        return eta

    def calc_amplitudes_with_dislocation_ring(self, eta=None,
                                                position=None,
                                                radius=None,
                                                normal_vector=[0, 0, 1],
                                                dislocation_type=1):
        """
        Inserts a  dislocation ring in the system corresponding to dislocation type

        """
        # TODO: Improve code documentation

        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 to insert a dislocation dipole.")

        if position == None:
            position = self.rmid

        if radius == None:
            radius = min(self.rmax)/3

        if eta == None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type - 1]

        theta = self.calc_angle_field_vortex_ring(
            radius=radius,
            position=position,
            normal_vector=normal_vector)

        for n in range(self.number_of_reciprocal_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * sn[n] * theta)

        return eta

    # PLOTTING FUNCTIONS


class PhaseFieldCrystal1DPeriodic(PhaseFieldCrystal):
    def __init__(self, nx, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal1DPeriodic'
        self.dim = 1

        # Default simulation parameters
        self.micro_resolution = [5]
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1]])
        self.q = np.array([[1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 1
        self.number_of_principal_reciprocal_modes = 1

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A = self.calc_initial_amplitudes()
        self.eta0 = [self.A]

        # Set the elastic constants
        self.el_mu = 3 * self.A ** 2
        self.el_lambda = 3 * self.A ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = (-3 * v * psi0 + np.sqrt(20 * t * v - 15 * v * r + 30 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) / (15 * v)
        return A

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2)


class PhaseFieldCrystal2DTriangular(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal2DTriangular'
        self.dim = 2
        self.nx = nx
        self.ny = ny

        # Default simulation parameters
        self.micro_resolution = [7, 12]
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]

        a0 = 2 * np.pi * 2 / np.sqrt(3)
        self.a = a0 * np.array([[1, 0], [1 / 2, np.sqrt(3) / 2], [1 / 2, -np.sqrt(3) / 2]])
        self.q = np.array([[np.sqrt(3) / 2, -1 / 2], [0, 1], [-np.sqrt(3) / 2, -1 / 2]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 3
        self.number_of_principal_reciprocal_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = np.sqrt(3) * a0 / self.micro_resolution[1]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A = self.calc_initial_amplitudes()
        self.eta0 = [self.A, self.A, self.A]

        # Set the elastic constants
        self.el_lambda = 3 * self.A ** 2
        self.el_mu = 3 * self.A ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = (-3 * v * psi0 + np.sqrt(20 * t * v - 15 * v * r + 30 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) / (15 * v)
        return A

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2)
    
    
    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array(
            [2*self.calc_gaussfilter_f()*np.fft.fftn(
            np.fft.ifftn((1-k2)*field_f)*np.fft.ifftn(self.dif[i]*k2*field_f)
            ) for i in range(self.dim)])



class PhaseFieldCrystal2DSquare(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal2DSquare'
        self.dim = 2
        self.nx = nx
        self.ny = ny

        # Default simulation parameters
        self.micro_resolution = [7, 7]
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0], [0, 1]])
        self.q = np.array([[1, 0], [0, 1], [1, -1], [1, 1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 4
        self.number_of_principal_reciprocal_modes = 2

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A, self.B = self.calc_initial_amplitudes()
        self.eta0 = [self.A, self.A, self.B, self.B]

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2
        self.el_mu = 16 * self.B ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        def equations(vars):
            A, B = vars
            eq1 = 12*psi0**2*A + 48*psi0*A*B + 36*A**3 + 72*A*B**2 + 4*A*r
            eq2 = 12*psi0**2*B + 24*psi0*A**2 + 36*B**3 + 72*A**2*B + 4*B*r
            return [eq1, eq2]


        A = 0
        B = 0

        for A0 in np.linspace(0,1,11):
            B0 = A0/2
            [A_tmp, B_tmp] = fsolve(equations, [A0, B0])

            if abs(A_tmp) > abs(A):
                A = A_tmp
                B = B_tmp

        return A, B

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2*(2-k2)**2)



class PhaseFieldCrystal3DBodyCenteredCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal3DBodyCenteredCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [7, 7, 7]
        self.psi0 = -0.325
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]
        self.zRes = nz * self.micro_resolution[2]

        a0 = 2 * np.pi * np.sqrt(2)
        self.a = a0 / 2 * np.array([[-1, 1, 1],
                                    [1, -1, 1],
                                    [1, 1, -1],
                                    [1, 1, 1]])

        self.q = np.array([[0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0]]) / np.sqrt(2)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 6
        self.number_of_principal_reciprocal_modes = 6

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A = self.calc_initial_amplitudes()
        self.eta0 = [self.A, self.A, self.A, self.A, self.A, self.A]

        # Set the elastic constants
        self.el_lambda = 4 * self.A ** 2
        self.el_mu = 4 * self.A ** 2
        self.el_gamma = - 4*self.A**2
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)


    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = -2/15 * psi0 + 1/15*np.sqrt(-5*r - 11*psi0**2)
        return A

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2)

class PhaseFieldCrystal3DFaceCenteredCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal3DBodyCenteredCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [11, 11, 11]
        self.psi0 = -0.325
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]
        self.zRes = nz * self.micro_resolution[2]

        a0 = 2 * np.pi * np.sqrt(3)
        self.a = a0 / 2 * np.array([[0, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0],
                                    [0, -1, 1],
                                    [-1, 0, 1],
                                    [-1, 1, 0]])

        self.q = np.array([[-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1],
                           [2, 0, 0],
                           [0, 2, 0],
                           [0, 0, 2]]) / np.sqrt(3)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 7
        self.number_of_principal_reciprocal_modes = 4

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A, self.B = self.calc_initial_amplitudes()
        self.eta0 = [self.A, self.A, self.A, self.A, self.B, self.B, self.B]

        # Set the elastic constants
        self.el_lambda = 32/81 * self.A ** 2
        self.el_mu = 32/81 * self.A ** 2
        self.el_gamma = 32/81 * (2*self.B**2 - self.A**2)
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        def equations(vars):
            A, B = vars
            eq1 = 27*A**2 + 36*B**2 + 18*B*psi0 + 3*psi0**2 + r
            eq2 = 72*A**2*(4*B + psi0) + 90*B**3 + 18*B*psi0**2 + 6*B*r
            return [eq1, eq2]

        A = 0
        B = 0

        for A0 in np.linspace(0, 1, 11):
            B0 = A0 / 2
            [A_tmp, B_tmp] = fsolve(equations, [A0, B0])

            if abs(A_tmp) > abs(A):
                A = A_tmp
                B = B_tmp

        return A, B

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2*(4/3-k2)**2)


class PhaseFieldCrystal3DSimpleCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal3DBodyCenteredCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [5, 5, 5]
        self.psi0 = -0.325
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1
        # TODO: Some of this code seems to be reprinted. Maybe it should be moved to the BaseSystem class? (Vidar 18.12.23)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]
        self.zRes = nz * self.micro_resolution[2]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]])

        self.q = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1],
                           [0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0],
                           [-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 13
        self.number_of_principal_reciprocal_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A, self.B, self.C = self.calc_initial_amplitudes()
        self.eta0 = [self.A, self.A, self.A,
                     self.B, self.B, self.B, self.B, self.B, self.B,
                     self.C, self.C, self.C, self.C]

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_mu = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_gamma = 32*self.A**2 - 16*self.B**2 - 256*self.C**2
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        def equations(vars):
            A, B, C = vars
            eq1 = 15*A**3 + 24*A**2*C + 24*B*C*(3*B + psi0) \
                    + 96*A*B**2 + 36*A*C**2 + 24*A*B*psi0 \
                    + 3*A*psi0**2 + A*r
            eq2 = 12*A*C*(6*B + psi0) + 6*A**2*(8*B + psi0) \
                    + 45*B**3 + 54*B*C**2 + 12*B**2*psi0 \
                    + 3*B*psi0**2 + B*r
            eq3 = 6*A**3 + 27*A**2*C + 18*A*B*(3*B + psi0) \
                    + C*(81*B**2 + 27*C**2 + 3*psi0**2 + r)

            return [eq1, eq2, eq3]

        A = 0
        B = 0
        C = 0

        for A0 in np.linspace(0, 0.2, 21):
            B0 = A0 / 2
            C0 = A0 / 4
            [A_tmp, B_tmp, C_tmp] = fsolve(equations, [A0, B0, C0])

            if abs(A_tmp) > abs(A):
                A = A_tmp
                B = B_tmp
                C = C_tmp

        return A, B, C

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2*(2-k2)**2*(3-k2)**2)