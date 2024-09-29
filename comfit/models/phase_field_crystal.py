import numpy as np
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint
import copy

# Plot functions
from comfit.plot.plot_field_plotly import plot_field_plotly
from comfit.plot.plot_field_matplotlib import plot_field_matplotlib

from comfit.tool.tool_print_in_color import tool_print_in_color


class PhaseFieldCrystal(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Initializes a system to simulate a Phase Field Crystal. This class is the base of the other
        phase field crystal models implemented in comfit

        Args:
            dimension: The dimension of the system.
            kwargs :keyword arguments to set additional parameters. See https://comfitlib.com/ClassPhaseFieldCrystal/

        Returns:
            The system object representing the phase field crystal simulation.
        """
        
        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.dislocation_charges = np.array(
            [[np.round(np.dot(an, qn) / (2 * np.pi), decimals=8) for qn in self.q] for an in self.a])

        # Elastic constant
        if ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)!=0:
            self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)
        else:
            self.el_nu = 0

        self.Phi = 2*sum(np.array(self.eta0)**2)

    def __str__(self):
        """
        Returns a string representation of the object.

        Args:
            None
        
        Returns: 
            A string representation of the object.
        """
        string = "-------------------------------\n"
        string += "Phase Field Crystal instance\n"
        string += "-------------------------------\n"
        string += "Type: " + self.type + "\n"
        if self.dim == 1:
            string += "(nx) = (" + str(self.nx) + ")\n"
        elif self.dim == 2:
            string += "(nx, ny) = (" + str(self.nx) + ", " + str(self.ny) + ")\n"
        elif self.dim == 3:
            string += "(nx, ny, nz) = (" + str(self.nx) + ", " + str(self.ny) + ", " + str(self.nz) + ")\n"
        string += "micro_resolution: " + str(self.micro_resolution) + "\n"
        string += "-------------------------------\n"
        string += "Parameters\n"
        string += "-------------------------------\n"
        string += "r: " + str(self.r) + "\n"
        string += "t: " + str(self.t) + "\n"
        string += "v: " + str(self.v) + "\n"
        string += "-------------------------------\n"
        string += "Proto amplitudes\n"
        string += "-------------------------------\n"
        string += "psi0: " + str(self.psi0) + "\n"
        if self.type == 'PhaseFieldCrystal1DPeriodic':
            string += "A: " + str(self.A) + "\n"
        elif self.type == 'PhaseFieldCrystal2DTriangular':
            string += "A: " + str(self.A) + "\n"
        elif self.type == 'PhaseFieldCrystal2DSquare':
            string += "A: " + str(self.A) + "\n"
            string += "B: " + str(self.B) + "\n"
        elif self.type == 'PhaseFieldCrystal3DBodyCenteredCubic':
            string += "A: " + str(self.A) + "\n"
        elif self.type == 'PhaseFieldCrystal3DFaceCenteredCubic':
            string += "A: " + str(self.A) + "\n"
            string += "B: " + str(self.B) + "\n"
        elif self.type == 'PhaseFieldCrystal3DSimpleCubic':
            string += "A: " + str(self.A) + "\n"
            string += "B: " + str(self.B) + "\n"
            string += "C: " + str(self.C) + "\n"
        string += "-------------------------------\n"

        return string
    
    #######################################################
    ############ CONFIGURATION FUNCTIONS #################
    #######################################################

    def conf_PFC_from_amplitudes(self, eta=None, rotation=None):
        """Configures the PFC from the amplitudes.

        Args:
            eta: The amplitudes to configure the PFC from.

        Returns:
            Configures self.psi and self.psi_f.
        """

        self.psi = self.calc_PFC_from_amplitudes(eta, rotation)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_advect_PFC(self,u):
        """Advects the PFC according to the displacement field u.

        Args:
            u: The displacement field to advect the PFC with.

        Returns:
            None, but updates the PFC.
        """

        self.psi = np.real(self.calc_advect_field(self.psi, u, self.psi_f))
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_apply_distortion(self, distortion):
        """Applies a distortion to the PFC.

        Args:
            distortion: The distortion to apply to the PFC.
        
        Returns:
            None, but updates the PFC.
        """

        self.bool_is_distorted = True

        distortion = np.array(distortion)

        if self.dim == 1:
            self.k[0] = self.k[0]/(1+distortion)
            self.dif[0] = self.dif[0]/(1+distortion)
            self.x = self.x*(1+distortion)
            self.xmax = self.xmax*(1+distortion)
            self.xmin = self.xmin*(1+distortion)
            self.size_x = self.size_x*(1+distortion)
            self.dx = self.dx*(1+distortion)

        elif self.dim == 2:
            # Creating 2D meshgrid
            X,Y = np.meshgrid(self.x.flatten(),self.y.flatten(), indexing='ij')

            # Strain matrix
            u_xx = distortion[0,0]
            u_xy = distortion[0,1]
            u_yx = distortion[1,0]
            u_yy = distortion[1,1]

            distortion_matrix = np.array([[1 + u_xx, u_xy    ],
                                          [u_yx,     1 + u_yy]])

            # Applying the strain
            X = distortion_matrix[0,0]*X + distortion_matrix[0,1]*Y
            Y = distortion_matrix[1,0]*X + distortion_matrix[1,1]*Y

            # Updating the x and y coordinates
            self.X = X
            self.Y = Y

            # Updating the k and dif vectors
            inverse_distortion_matrix = np.linalg.inv(distortion_matrix)

            self.k[0] = self.k[0]*inverse_distortion_matrix[0,0] + self.k[1]*inverse_distortion_matrix[0,1]
            self.k[1] = self.k[0]*inverse_distortion_matrix[1,0] + self.k[1]*inverse_distortion_matrix[1,1]

            self.dif[0] = 1j*self.k[0]
            self.dif[1] = 1j*self.k[1]

            # Updating the dx and dy
            original_dx = self.dx
            original_dy = self.dy
            self.dx = np.sqrt((distortion_matrix[0,0]*original_dx)**2 + (distortion_matrix[0,1]*original_dy)**2)
            self.dy = np.sqrt((distortion_matrix[1,0]*original_dx)**2 + (distortion_matrix[1,1]*original_dy)**2)

            # Updating the size_x and size_y
            original_size_x = self.size_x
            original_size_y = self.size_y
            self.size_x = np.sqrt((distortion_matrix[0,0]*original_size_x)**2 + (distortion_matrix[0,1]*original_size_y)**2)
            self.size_y = np.sqrt((distortion_matrix[1,0]*original_size_x)**2 + (distortion_matrix[1,1]*original_size_y)**2)

            # Updating the x and y coordinate limits
            self.xmax = np.max(X)
            self.xmin = np.min(X)

            self.ymax = np.max(Y)
            self.ymin = np.min(Y)

            self.volume = abs(np.linalg.det(distortion_matrix))*self.volume


        else:
            raise ImplementationError("Applied distortion is not yet implemented for 3 dimensions.")

    def conf_strain_to_equilibrium(self):
        """Configures 
        Strains the pfc to equilibrium by adjusting the position variables and k-space variables.

        Args:
            None

        Returns:
            None, configures the PFC.
        """
        
        # Evolve PFC to reach initial state of equilibrium (100 time steps ok)
        number_of_initial_steps = 100
        self.evolve_PFC(number_of_initial_steps, suppress_output=True)
        number_of_steps = 10 # 10 steps is enough to reach equilibrium

        # Calculate average free energy to compare
        average_free_energy = self.calc_free_energy()/self.volume
        # print('Average free energy at initial state: ', average_free_energy)

        strain = 0
        strain_increment = 0.00001
        strain += strain_increment

        # Unaltered k-vectors
        k0 = self.k[0].copy()
        dx0 = self.dx
        x0 = self.x.copy()
        xmin0 = self.xmin
        xmax0 = self.xmax
        size_x0 = self.size_x
        a00 = self.a0

        if self.dim > 1:
            k1 = self.k[1].copy()
            dy0 = self.dy
            y0 = self.y.copy()
            ymin0 = self.ymin
            ymax0 = self.ymax
            size_y0 = self.size_y

            
        if self.dim > 2:
            k2 = self.k[2].copy()
            dz0 = self.dz
            z0 = self.z.copy()
            zmin0 = self.zmin
            zmax0 = self.zmax
            size_z0 = self.size_z

        volume0 = self.volume
        dV0 = self.dV
        
        def update_lengths(self, strain):
            self.k[0] = k0/(1+strain)
            self.dx = dx0*(1+strain)
            self.x = x0*(1+strain)
            self.xmax = xmax0*(1+strain)
            self.xmin = xmin0*(1+strain)
            self.size_x = size_x0*(1+strain)
            self.a0 = a00*(1+strain)
            
            if self.dim > 1:
                self.k[1] = k1/(1+strain)
                self.dy = dy0*(1+strain)
                self.y = y0*(1+strain)
                self.ymax = ymax0*(1+strain)
                self.ymin = ymin0*(1+strain)
                self.size_y = size_y0*(1+strain)
                
            if self.dim > 2:
                self.k[2] = k2/(1+strain)
                self.dz = dz0*(1+strain)
                self.z = z0*(1+strain)
                self.zmax = zmax0*(1+strain)
                self.zmin = zmin0*(1+strain)
                self.size_z = size_z0*(1+strain)

            self.dV = self.dx*self.dy*self.dz
            self.volume = self.size_x*self.size_y*self.size_z

        update_lengths(self,strain)

        # Evolve PFC
        self.evolve_PFC(number_of_steps, suppress_output=True)

        # Calculate free energy to compare
        average_free_energy_tmp = self.calc_free_energy()/self.volume

        # print('Free energy at strain: ', strain, ' is: ', free_energy_tmp)

        positive_strain_required = False
        while average_free_energy_tmp < average_free_energy:

            positive_strain_required = True

            average_free_energy = average_free_energy_tmp

            strain += strain_increment

            update_lengths(self,strain)
            
            self.evolve_PFC(number_of_steps, suppress_output=True)

            average_free_energy_tmp = self.calc_free_energy()/self.volume
            # print('Average free energy at strain: ', strain, ' is: ', average_free_energy_tmp)
        
        if positive_strain_required:
            # Going one back to get the lowest free energy
            final_strain = strain - strain_increment
            update_lengths(self,strain)
            # tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')

        else: #negative strain required

            strain = - strain_increment

            update_lengths(self,strain)
            self.evolve_PFC(number_of_steps, suppress_output=True)
            average_free_energy_tmp = self.calc_free_energy()/self.volume
            print('Average Free energy at strain: ', strain, ' is: ', average_free_energy_tmp)

            while average_free_energy_tmp < average_free_energy:

                average_free_energy = average_free_energy_tmp

                strain -= strain_increment

                update_lengths(self, strain)

                self.evolve_PFC(number_of_steps, suppress_output=True)

                average_free_energy_tmp = self.calc_free_energy()/self.volume

                # print('Average free energy at strain: ', strain, ' is: ', average_free_energy_tmp)

            # Going one back to get the lowest free energy
            final_strain = strain + strain_increment
            update_lengths(self, final_strain)
            # tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')   

        return final_strain

    def conf_create_polycrystal(self, type, **kwargs):
        """Creates a polycrystal.

        Args:
            type: The type of polycrystal to create.
            kwargs: Additional arguments for the polycrystal creation. Including:
                relaxation_time: The relaxation time to use for the polycrystal creation.
        
        Returns:
            None, but updates the PFC.
        """
        # First type of polycrystal, see documentation.
        if type == 'circular':
            # This creates a standard orientation of the crystal
            self.conf_PFC_from_amplitudes()

            rotation = kwargs.get('rotation', np.pi/6)
            position = kwargs.get('position', self.rmid)
            radius = kwargs.get('radius', self.size_min/4)

            # Create the rotated field
            psi_rotated = self.calc_PFC_from_amplitudes(rotation=[0,0,rotation])

            # Set the rotated field in the inclusion region
            region  = self.calc_region_disk(position, radius)
            self.psi[region] = psi_rotated[region]
            self.psi_f = sp.fft.fftn(self.psi)

            # Smooth the interface
            relaxation_time = kwargs.get('relaxation_time', 10)
            self.evolve_PFC(round(relaxation_time/self.dt))

        elif type == 'four_grain':
            if self.dim == 1:
                raise Exception("Polycrystal type four_grain is not valid for 1 dimension.") 
                
            self.psi = self.calc_PFC_from_amplitudes(self.eta0)
            
            l1  = self.y>1/6*self.ymax+(4/6*self.ymax)/(1/1*self.xmax)*self.x
            l2  = self.y>1/2*self.ymax+(1/2*self.ymax)/(2/3*self.xmax)*self.x
            l3  = self.y>1/2*self.ymax-(1/2*self.ymax)/(2/3*self.xmax)*self.x
            l4  = self.y>5/6*self.ymax-(4/6*self.ymax)/(1/1*self.xmax)*self.x
            l5  = self.x>self.xmax/4
            l6  = self.y>5/4*self.ymax-(3/4*self.ymax)/(1/1*self.xmax)*self.x
            l7  = self.x>self.xmax/2
            l8  = self.y<1/2*self.ymax+(1/2*self.ymax)/(2/3*self.xmax)*self.x
            l9  = self.x>3/4*self.xmax
            l10 = self.y<-1/4*self.ymax+(3/4*self.ymax)/(1/1*self.xmax)*self.x

            zdir = 1 if self.dim == 2 else np.ones((1,1,self.zRes))

            pfcRotated = self.calc_PFC_from_amplitudes(self.eta0, rotation=[0,0,-22.5/180*np.pi])
            region = np.bool_((l4*~(l7)*~(l8) + ~(l1)*~(l3)*~(l7))*zdir)
            self.psi[region] = pfcRotated[region]

            pfcRotated = self.calc_PFC_from_amplitudes(self.eta0, rotation=[0,0,22.5/180*np.pi])
            region = np.bool_((l1*l6*l7 + l7*l10*~(l4))*zdir)
            self.psi[region] = pfcRotated[region]

            pfcRotated = self.calc_PFC_from_amplitudes(self.eta0, rotation=[0,0,45/180*np.pi])
            region = np.bool_((l1*~(l4)*~(l5) + ~(l1)*l4*l9)*zdir)
            self.psi[region] = pfcRotated[region]

            self.psi_f = sp.fft.fftn(self.psi)

            relaxation_time = kwargs.get('relaxation_time', 10)
            
            self.evolve_PFC(round(relaxation_time/self.dt))


    #######################################################
    #######################################################
    ############### EVOLUTION FUNCTIONS ###################
    #######################################################
    #######################################################


    #######################################################
    #################### STANDARD #########################
    #######################################################
    def evolve_PFC(self, number_of_steps, method='ETD2RK', suppress_output=False):
        """Evolves the PFC according to classical PFC dynamics.

        Args:
            number_of_steps: The number of steps to evolve the PFC.
            method: The method to use for the evolution (default: ETD2RK).
        
        Returns:
            Updates self.psi and self.psi_f
        """

        if self.type_of_evolution == 'conserved':
            omega_f = -self.calc_k2()*(self.r + self.calc_L_f()**2)
            non_linear_evolution_function_f = self.calc_nonlinear_evolution_function_conserved_f
        elif self.type_of_evolution == 'unconserved':
            omega_f = -(self.r + self.calc_L_f()**2)
            non_linear_evolution_function_f = self.calc_nonlinear_evolution_function_unconserved_f

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the PFC (conserved)', disable=suppress_output):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_conserved_f,
                                          self.psi, self.psi_f)

            # These steps seem to be necessary for numerical stability (Vidar 18.12.23)
            self.psi = np.real(self.psi)
            self.psi_f = sp.fft.fftn(self.psi)

    # Nonlinear part conserved
    def calc_nonlinear_evolution_function_conserved_f(self, psi, t):
        return -self.calc_k2()*sp.fft.fftn(self.t * psi ** 2 + self.v * psi ** 3)

    # Non-linear part unconserved
    def calc_nonlinear_evolution_function_unconserved_f(self, psi, t):
        return -sp.fft.fftn(self.t * psi ** 2 + self.v * psi ** 3)
    #######################################################

    #######################################################
    ##### MECHANICAL EQUILIBRIUM (CONSERVED) ##############
    #######################################################
    def evolve_PFC_mechanical_equilibrium(self, time, Delta_t = 10, method='ETD2RK'):
        """Evolves the PFC in mechanical equilibrium. 

        Args:
            number_of_steps: The number of steps to evolve the PFC
            Delta_t: The time step for the mechanical equilibrium evolution
            method: The method to use for the evolution

        Returns:
            Updates self.psi and self.psi_f
        """
        
        number_of_steps = round(time/self.dt)
        number_of_steps_per_iteration = round(Delta_t/self.dt)
        number_of_iterations = round(number_of_steps/number_of_steps_per_iteration)

        for n in tqdm(range(number_of_iterations), desc='Evolving the PFC in mechanical equilibrium'):
            self.conf_advect_PFC(self.calc_displacement_field_to_equilibrium())
            self.evolve_PFC(number_of_steps_per_iteration, method)
    #######################################################

    #######################################################
    ######### HYDRODYNAMIC (CONSERVED) ####################
    #######################################################
    def evolve_PFC_hydrodynamic(self, number_of_steps, 
                                method = 'ETD2RK',
                                gamma_S = 2**-4,
                                rho0 = 2**-4):
        """Evolves the PFC according to hydrodynamic PFC dynamics.

        This requires introducing a velocity field. 
        If psi does not contain this field, it is added to the components psi[1],psi[2],psi[3].
        
        Args:
            number_of_steps: The number of steps to evolve the PFC.
            method: The method to use for the evolution (default: ETD2RK).
            gamma_S: The surface tension coefficient.
            rho0: The mass density.

        Returns:
            Updates self.psi and self.psi_f
        """
        
        if hasattr(self,'bool_has_velocity_field'):
            pass
        else:
            self.bool_has_velocity_field = True
            self.psi = np.array([self.psi]+[np.zeros_like(self.psi)]*self.dim)
            self.psi_f = np.array([self.psi_f]+[np.zeros_like(self.psi_f)]*self.dim)
            # print("psi shape", self.psi.shape)

            if not hasattr(self,'external_force_density_f'):
                self.external_force_density_f = np.zeros([self.dim] + self.dims, dtype=complex)

        self.gamma_S = gamma_S
        self.rho0 = rho0

        omega_f = self.calc_omega_hydrodynamic_f()
        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the hydrodynamic PFC'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_hydrodynamic_evolution_function_f,
                                          self.psi, self.psi_f)
            
            self.psi = np.real(self.psi)
            self.psi_f = sp.fft.fftn(self.psi, axes = (range ( - self.dim , 0) ) )

    # Linear part
    def calc_omega_hydrodynamic_f(self):
        """Calculates the hydrodynamic evolution function omega_f.

        Args:
            None

        Returns:
            The hydrodynamic evolution function omega_f.
        """
        k2 = self.calc_k2()
        return np.array([self.calc_omega_f()]+[-self.gamma_S/self.rho0*k2]*self.dim)

    # Nonlinear part
    def calc_nonlinear_hydrodynamic_evolution_function_f(self, field, t):
        """Calculates the hydrodynamic evolution function of the PFC.

        Args:
            field: The field to calculate the evolution function of.
            t: The time.

        Returns:
            The nonlinear evolution function for the hydrodynamic PFC.
        """

        field_f = sp.fft.fftn(field, axes =( range ( - self . dim , 0) ))

        k2 = self.calc_k2()

        N0_f = -k2*sp.fft.fftn(self.t * field[0] ** 2 + self.v * field[0] ** 3) \
            - sp.fft.fftn(sum([field[i+1]*sp.fft.ifftn(self.dif[i]*field_f[0]) for i in range(self.dim)]))
        
        force_density_f = self.calc_stress_divergence_f(field_f[0])

        return np.array([N0_f] + [1/self.rho0*(force_density_f[i]+self.external_force_density_f[i]) for i in range(self.dim)])
    #######################################################

    #######################################################
    #######################################################
    ############## CALCULATION FUNCTIONS ##################
    #######################################################
    #######################################################
    
    def calc_strained_amplitudes(self):
        """ Straines the PFC to equilibrium and returns the amplitudes.
    
        Args:
            None

        Returns:
            The amplitudes of the strained PFC.
        """
        tool_print_in_color('Proto amplitudes and elastic constants', 'blue')
        tool_print_in_color('---', 'blue')
        print(f'Proto psi0: {self.psi0:.02f}')
        print(f'Proto A: {self.A:.05f}')
        if self.type in ['PhaseFieldCrystal2DSquare','PhaseFieldCrystal3DFaceCenteredCubic','PhaseFieldCrystal3DSimpleCubic']:
            print(f'Proto B: {self.B:.05f}')
        if self.type in ['PhaseFieldCrystal3DSimpleCubic']:
            print(f'Proto C: {self.C:.05f}')

        print(f'Proto mu: {self.el_mu:.05f}')
        print(f'Proto lambda: {self.el_lambda:.05f}')
        print(f'Proto gamma: {self.el_gamma:.05f}')

        # Strain PFC to equilibrium
        self.conf_PFC_from_amplitudes()
        final_strain = self.conf_strain_to_equilibrium()
        tool_print_in_color('Amplitudes after strain', 'blue')
        tool_print_in_color('---', 'blue')
        print(f'Equilibrium strain: {final_strain:.05f}')
        print(f'Equilibrium q-vector: {1/(1+final_strain):.05f}')


        psi0 = np.mean(self.psi)
        print(f'Eq. psi0: {psi0:.02f}')

        eta = self.calc_demodulate_PFC()
        A = np.mean(np.real(eta[0]))
        number_of_independent_amplitudes = 1

        if self.type == 'PhaseFieldCrystal2DSquare':
            number_of_independent_amplitudes = 2
            B = np.mean(np.real(eta[2]))
        
        if self.type == 'PhaseFieldCrystal3DFaceCenteredCubic':
            number_of_independent_amplitudes = 2
            B = np.mean(np.real(eta[4]))

        if self.type == 'PhaseFieldCrystal3DSimpleCubic':
            number_of_independent_amplitudes = 3
            B = np.mean(np.real(eta[3]))
            C = np.mean(np.real(eta[9]))

        if number_of_independent_amplitudes >= 1:
            print(f'Eq. A: {A:.05f}')
            print('A eq./A proto: {:.05f}'.format(A/self.A))
        if number_of_independent_amplitudes >= 2:
            print(f'Eq. B: {B:.05f}')
            print('Ratio B eq./B proto: {:.05f}'.format(B/self.B))
        if number_of_independent_amplitudes >= 3:
            print(f'C: {C:.05f}')
            print('Ratio C eq./C proto: {:.05f}'.format(C/self.C))
        

        #Finding elastic constants, starting with mu
        pfc_strained = copy.deepcopy(self)
        f0 = pfc_strained.calc_free_energy()/pfc_strained.volume

        shear_strain=-0.001
        pfc_strained.conf_apply_distortion(np.array([[0,shear_strain],[0.0,0.0]]))
        f = pfc_strained.calc_free_energy()/pfc_strained.volume

        # TO be refined
        mu = 2*(f-f0)/(shear_strain**2)
        print(f'Eq. mu: {mu:.05f}')
        print('Ratio mu eq./mu proto: {:.05f}'.format(mu/self.el_mu))


        # gamma
        pfc_strained = copy.deepcopy(self)
        f0 = pfc_strained.calc_free_energy()/pfc_strained.volume

        compression_strain1=0.0001
        pfc_strained.conf_apply_distortion(np.array([[compression_strain1,0],[0,-compression_strain1]]))
        # print('Volume ratio after compression: {:.05f}'.format(pfc_strained.volume/self.volume))
        f1 = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        pfc_strained = copy.deepcopy(self)
        compression_strain2=0.0002
        pfc_strained.conf_apply_distortion(np.array([[compression_strain2,0],[0,-compression_strain2]]))
        f2 = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        # gamma = (f-f0 - 2*mu*compression_strain**2)/(compression_strain**2)
        gamma = (f1 - f2)/(compression_strain1**2-compression_strain2**2) - 2*mu
        print(f'Eq. gamma: {gamma:.05f}')
        print('Ratio gamma eq./gamma proto: {:.05f}'.format(gamma/self.el_gamma))


    def calc_PFC_free_energy_density_and_chemical_potential(self,field=None,field_f=None):
        """Calculates the free energy density and chemical potential of the PFC.

        Args:
            field: The field to calculate the free energy density and chemical potential of.
            field_f: The Fourier transform of the field.

        Returns:
            The free energy density and chemical potential of the PFC.
        """

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

        Lpsi = sp.fft.ifftn(self.calc_L_f()*psi_f)

        free_energy_density = 1/2*Lpsi**2 \
            + 1/2*self.r*psi2 + 1/3*self.t*psi3 + 1/4*self.v*psi4
        
        # print("Lpsi shape",self.calc_Lpsi(psi_f).shape)

        L2psi = sp.fft.ifftn(self.calc_L_f()**2*psi_f)

        chemical_potential =L2psi*psi \
            + self.r*psi + self.t*psi2 + self.v*psi3
        
        # print("L2psi shape",self.calc_L2psi(psi_f).shape)
        
        # print("Free energy density shape", free_energy_density.shape)
        # print("Chem pot shape", chemical_potential.shape)
        return np.real(free_energy_density), np.real(chemical_potential)
    

    def calc_displacement_field_to_equilibrium(self):
        """Calculates the displacement field needed to put the PFC in mechanical equilibrium.

        Args:
            None
        
        Returns:
            The displacement field u
        """

        # Calculate the stress divergence
        g_f = self.calc_stress_divergence_f()

        # Calculate the displacement field
        u_f = np.zeros([self.dim] + self.dims, dtype = complex)

        k2 = self.calc_k2()
        k = np.sqrt(k2)

        kappa = np.zeros([self.dim] + self.dims)
        for l in range(self.dim):
            kappa[l] = self.k[l]/k

        second_term_factor = (self.el_mu + self.el_lambda)/(1+sum([(self.el_mu + self.el_lambda)/(self.el_mu + self.el_gamma*kappa[l]**2)*kappa[l]**2 for l in range(self.dim)]))

        # print(second_term_factor)

        for i in range(self.dim):
            for j in range(self.dim):
                #Calculate the Green's function
                Gij_f = int(i==j)/(self.el_mu + self.el_gamma*kappa[i]**2) \
                - kappa[i]*kappa[j]/((self.el_mu + self.el_gamma*kappa[i]**2)*(self.el_mu + self.el_gamma*kappa[j]**2))*second_term_factor

                u_f[i] += 1/k2*Gij_f*g_f[j]
            
            # Set the zero mode to zero
            u_f[i][self.zero_index] = 0

            # print(u_f[i])

        return np.real(sp.fft.ifftn(u_f, axes = (range ( - self . dim , 0) )))

        
    # Initial configuration methods
    def calc_amplitudes_with_dislocation(self, eta=None, x=None, y=None, dislocation_type=1):
        """ Calculates the amplitudes with a single point dislocation inserted.

        Args:
            eta: The amplitudes to insert the dislocation in.
            x: The x-coordinate of the dislocation.
            y: The y-coordinate of the dislocation.
            dislocation_type: The dislocation type to insert.
        
        Returns:
            The amplitudes containing the dislocation
        """

        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a single point dislocation.")

        if x == None:
            x = self.xmid
        if y == None:
            y = self.ymid
        if eta == None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type - 1]
        for n in range(self.number_of_reciprocal_lattice_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * self.calc_angle_field_single_vortex([x, y], charge=sn[n]))

        return eta

    def calc_amplitudes_with_dislocation_dipole(self, eta=None,
                                                x1=None, y1=None,
                                                x2=None, y2=None,
                                                dislocation_type=1):
        """
        Inserts a  dislocation dipole in the system corresponding to dislocation type and its negative

        Args:
            eta: The amplitudes to insert the dislocation dipole in
            x1: The x-coordinate of the first dislocation
            y1: The y-coordinate of the first dislocation
            x2: The x-coordinate of the second dislocation
            y2: The y-coordinate of the second dislocation
            dislocation_type: The dislocation type to insert
        
        Returns:
            The amplitudes with the dislocation dipole inserted
        """

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
            eta = self.eta0.reshape(self.number_of_reciprocal_lattice_modes,1,1)\
                *np.ones([self.number_of_reciprocal_lattice_modes] + self.dims, dtype=complex)

        sn = self.dislocation_charges[dislocation_type - 1]

        theta = self.calc_angle_field_vortex_dipole(
            dipole_vector=[x2 - x1, y2 - y1],
            dipole_position=[(x1 + x2) / 2, (y1 + y2) / 2])

        for n in range(self.number_of_reciprocal_lattice_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * sn[n] * theta)

        return eta

    def calc_amplitudes_with_dislocation_ring(self, eta=None,
                                                position=None,
                                                radius=None,
                                                normal_vector=[0, 0, 1],
                                                dislocation_type=1):
        """
        Inserts a  dislocation ring in the system corresponding to dislocation type.

        Args:
            eta: The amplitudes to insert the dislocation ring in
            position: The position of the dislocation ring
            radius: The radius of the dislocation ring
            normal_vector: The normal vector of the dislocation ring
            dislocation_type: The dislocation type to insert

        Returns:
            The amplitudes with the dislocation ring inserted

        """

        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 to insert a dislocation dipole.")

        if position == None:
            position = self.rmid

        if radius == None:
            radius = min(self.rmax)/3

        if eta == None:
            eta = self.eta0.reshape(self.number_of_reciprocal_lattice_modes,1,1,1)\
                *np.ones([self.number_of_reciprocal_lattice_modes] + self.dims, dtype=complex)

        sn = self.dislocation_charges[dislocation_type - 1]

        theta = self.calc_angle_field_vortex_ring(
            radius=radius,
            position=position,
            normal_vector=normal_vector)

        for n in range(self.number_of_reciprocal_lattice_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j * sn[n] * theta)

        return eta
    
    def calc_PFC_from_amplitudes(self, eta=None, rotation=None):
        """Calculates the PFC from the amplitudes.

        Args:
            eta: The amplitudes to calculate the PFC from.
            rotation : The rotation of the PFC

        Returns:
            The PFC.
        """ 

        psi = self.psi0

        if rotation is None:
            rotation = [0,0,0]
        else:
            if isinstance(rotation, float):
                rotation = [0,0,rotation]

        rotation = sp.spatial.transform.Rotation.from_rotvec(rotation)
        rotation_matrix = rotation.as_matrix()

        # Rotate q-vectors to the new orientation
        if self.dim == 1:
            q = self.q # No rotation in 1D
        if self.dim == 2:
            q = (rotation_matrix[0:2,0:2]@self.q.transpose()).transpose()
        elif self.dim == 3:
            q = (rotation_matrix[0:3,0:3]@self.q.transpose()).transpose()

        if eta is None:
            eta = self.eta0

        for n in range(self.number_of_reciprocal_lattice_modes):

            if self.dim == 1:
                psi += 2 * eta[n] * np.exp(1j * q[n][0] * self.x)

            elif self.dim == 2:
                psi += 2 * eta[n] * np.exp(1j * (q[n][0] * self.x + q[n][1] * self.y))

            elif self.dim == 3:
                psi += 2 * eta[n] * np.exp(1j * (q[n][0] * self.x + q[n][1] * self.y + q[n][2] * self.z))

        return np.real(psi)

    def calc_demodulate_PFC(self):
        """Demodulates the PFC.

        Args:
            None

        Returns:
            The amplitudes corresponding to the demodulated PFC.
        """

        eta = np.zeros([self.number_of_primary_reciprocal_lattice_modes] + self.dims, 
                       dtype=complex)

        Gaussian_filter_f = self.calc_Gaussian_filter_f()

        order_parameter = self.psi if self.psi.ndim == self.dim else self.psi[0]

        if self.dim == 2:
                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(order_parameter*np.exp(
                        -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y)))

        elif self.dim == 3:
            for n in range(self.number_of_primary_reciprocal_lattice_modes):
                eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(order_parameter*np.exp(
                    -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y - 1j*self.q[n][2]*self.z  
                    )))
                
        return eta

    def calc_stress_tensor(self):
        """Calculates the stress of the phase-field crystal.

        Args:
            None

        Returns:
            The stress tensor
        """
        
        stress = self.calc_stress_tensor_microscopic()
        Gaussian_filter_f = self.calc_Gaussian_filter_f()

        # Number of independent stress components
        if self.dim == 1:
            number_of_stress_components = 1
        elif self.dim == 2:
            number_of_stress_components = 3
        elif self.dim == 3:
            number_of_stress_components = 6

        # Coarse-grain the microscopic stress
        for n in range(number_of_stress_components):
            stress[n] = Gaussian_filter_f*sp.fft.fftn(stress[n])

        return np.real(sp.fft.ifftn(stress, axes = (range( - self.dim , 0) )))

    def calc_structure_tensor(self):
        """Calculates the structure tensor of the phase-field crystal.

        Args:
            None

        Returns:
            The structure tensor
        """
        # Calculate the gradient
        diPsi = np.zeros([self.dim] + self.dims, dtype=complex)
        for i in range(self.dim):
            diPsi[i] = self.dif[i]*self.psi_f
        diPsi = np.real(sp.fft.ifftn(diPsi, axes = (range( - self.dim , 0) )))

        if self.dim == 1:
            number_of_independent_strain_components = 1
        elif self.dim == 2:
            number_of_independent_strain_components = 3
        elif self.dim == 3:
            number_of_independent_strain_components = 6

        structure_tensor = np.zeros([number_of_independent_strain_components] + self.dims, dtype=complex)

        Gaussian_filter_f = self.calc_Gaussian_filter_f()
        if self.dim == 1:
            structure_tensor[0] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[0])
        elif self.dim == 2:
            structure_tensor[0] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[0])
            structure_tensor[1] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[1])
            structure_tensor[2] = Gaussian_filter_f*sp.fft.fftn(diPsi[1]*diPsi[1])
        elif self.dim == 3:
            structure_tensor[0] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[0])
            structure_tensor[1] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[1])
            structure_tensor[2] = Gaussian_filter_f*sp.fft.fftn(diPsi[0]*diPsi[2])
            structure_tensor[3] = Gaussian_filter_f*sp.fft.fftn(diPsi[1]*diPsi[1])
            structure_tensor[4] = Gaussian_filter_f*sp.fft.fftn(diPsi[1]*diPsi[2])
            structure_tensor[5] = Gaussian_filter_f*sp.fft.fftn(diPsi[2]*diPsi[2])
        
        structure_tensor = np.real(sp.fft.ifftn(structure_tensor, axes = (range( - self.dim , 0) )))

        return structure_tensor
            


    def calc_strain_tensor(self):
        """Calculates the strain of the phase-field crystal.

        Args:
            None

        Returns:
            The strain tensor
        """
        
        strain = -self.dim/(2*self.Phi)*self.calc_structure_tensor()

        if self.dim == 1:
            strain = 1/2 + strain
        elif self.dim == 2:
            strain[0] = 1/2 + strain[0]
            strain[2] = 1/2 + strain[2]
        elif self.dim == 3:
            strain[0] = 1/2 + strain[0]
            strain[3] = 1/2 + strain[3]
            strain[5] = 1/2 + strain[5]
        
        return strain
        



    def calc_dislocation_density(self, eta = None):
        """Calculates the dislocation density

        Args:
            eta: The amplitudes to calculate the dislocation density from.

        Returns:
            The dislocation density
        """

        if eta is None:
            eta = self.calc_demodulate_PFC()

        if self.dim == 2:
            alpha = np.zeros([2] + self.dims)

            for n in range(self.number_of_primary_reciprocal_lattice_modes):
                D = self.calc_determinant_field([np.real(eta[n]), np.imag(eta[n])])
                alpha[0] += D*self.q[n,0]
                alpha[1] += D*self.q[n,1]

        elif self.dim == 3:
            alpha = np.zeros([3,3] + self.dims)

            for n in range(self.number_of_primary_reciprocal_lattice_modes):
                D = self.calc_determinant_field([np.real(eta[n]), np.imag(eta[n])])
                for i in range(3):
                    alpha[i][0] += D[i]*self.q[n,0]
                    alpha[i][1] += D[i]*self.q[n,1]
                    alpha[i][2] += D[i]*self.q[n,2]

        alpha = 2*self.dim/(self.number_of_primary_reciprocal_lattice_modes*self.A**2)*alpha

        return alpha
        

    def calc_dislocation_nodes(self):
        """ Calculates the dislocation nodes

        Args:
            None
        
        Returns:
            The dislocation nodes
        """

        # At the moment (24.09.24 - Vidar), this function only works for unstrained PFCs
        if hasattr(self,'bool_is_distorted') and self.bool_is_distorted:
            raise Exception("Dislocation nodes cannot be calculated for distorted PFCs.")

        alpha = self.calc_dislocation_density()

        if self.dim == 2:
            # self.plot_field(np.sqrt(alpha[0]**2 + alpha[1]**2))
            # plt.show()
            dislocation_nodes = self.calc_defect_nodes(np.sqrt(alpha[0]**2 + alpha[1]**2),
                                                       charge_tolerance=0.2*self.a0)
            
            for dislocation_node in dislocation_nodes:
                Burgers_vector = np.array([
                    alpha[0][dislocation_node['position_index']], 
                    alpha[1][dislocation_node['position_index']]
                    ])
                
                # print("Burgers vector:", Burgers_vector)
                
                # Find the Burgers vector 
                biggest_overlap = 0
                for a in self.a:
                    for sign in [-1, 1]:
                        overlap = np.dot(Burgers_vector, sign*a)
                        # print("Burgers vector",Burgers_vector)
                        # print("a-vector",a*sign)
                        # print("overlap",overlap)
                        # print("biggest overlap", biggest_overlap)
                        
                        if overlap > biggest_overlap:
                            # print("tjobing")
                            biggest_overlap = overlap
                            dislocation_node['Burgers_vector'] = sign*a
                
        elif self.dim == 3:
            dislocation_nodes = self.calc_defect_nodes(
                np.sqrt(
                      alpha[0][0]**2 + alpha[0][1]**2 + alpha[0][2]**2 \
                    + alpha[1][0]**2 + alpha[1][1]**2 + alpha[1][2]**2 \
                    + alpha[2][0]**2 + alpha[2][1]**2 + alpha[2][2]**2
                )
            )

            for dislocation_node in dislocation_nodes:
                alpha_tensor = np.array([
                    [alpha[0][0][dislocation_node['position_index']], alpha[0][1][dislocation_node['position_index']], alpha[0][2][dislocation_node['position_index']]],
                    [alpha[1][0][dislocation_node['position_index']], alpha[1][1][dislocation_node['position_index']], alpha[1][2][dislocation_node['position_index']]],
                    [alpha[2][0][dislocation_node['position_index']], alpha[2][1][dislocation_node['position_index']], alpha[2][2][dislocation_node['position_index']]]
                ])
                # print("alpha-tensor:",alpha_tensor)

                U, S, V = np.linalg.svd(alpha_tensor)
                tangent_vector = U[:,0]
                Burgers_vector = V[0,:]

                # print("U:",U)
                # print("S:",S)   
                # print("V:",V)

                # print("tangent vector:", tangent_vector)
                # print("Burgers vector:", Burgers_vector)

                # Find the Burgers vector 
                biggest_overlap = 0
                tangent_vector_sign = np.nan
                # print("Finding new Burgers vector...")
                for a in self.a:
                    for sign in [-1, 1]:
                        # print("Burgers vector",Burgers_vector)
                        # print("a-vector",a*sign)
                        # print("overlap",np.dot(Burgers_vector,sign*a))
                        overlap = np.dot(Burgers_vector, sign*a)
                        if overlap > biggest_overlap:
                            biggest_overlap = overlap
                            # print("New biggest overlap", overlap)
                            # print("Setting new BUrgers vector",a)
                            dislocation_node['Burgers_vector'] = a
                            tangent_vector_sign = sign

                dislocation_node['tangent_vector'] = tangent_vector_sign*tangent_vector

        return dislocation_nodes

    def calc_orientation_field(self):
        """Calculates the orientation field of the phase-field crystal.

        Args:
            None
        
        Returns:
            An orientation field, which is a vector field specifying the orientation of the crystal.
        """

        order_parameter = self.psi if self.psi.ndim == self.dim else self.psi[0]

        eta = np.zeros([self.number_of_primary_reciprocal_lattice_modes] + self.dims, 
                       dtype=complex)
        Gaussian_filter_f = self.calc_Gaussian_filter_f()

        if self.dim == 2:
            if self.type == 'PhaseFieldCrystal2DTriangular':
                wrap_factor = 3
            elif self.type == 'PhaseFieldCrystal2DSquare':
                wrap_factor = 2

            resolution = 16
            orientation_field = np.zeros((2,self.xRes,self.yRes))
            angles = np.arange(resolution)/resolution*np.pi/wrap_factor

            for angle in angles:
                rotation = sp.spatial.transform.Rotation.from_rotvec([0,0,angle])
                rotation_matrix = rotation.as_matrix()
                q = (rotation_matrix[0:2,0:2]@self.q.transpose()).transpose()
                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(order_parameter*np.exp(
                            -1j*q[n][0]*self.x - 1j*q[n][1]*self.y)))
                
                Phi2 = np.zeros(self.dims)
                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    Phi2 += abs(eta[n])**2
                
                orientation_field[0] += Phi2*np.cos(2*wrap_factor*angle)
                orientation_field[1] += Phi2*np.sin(2*wrap_factor*angle)

            return orientation_field
        
        elif self.dim == 3:

            # 1153 points
            xRes = 17
            yRes = 17
            zRes = 9

            x = np.linspace(-np.pi,np.pi, xRes)
            y = np.linspace(-np.pi,np.pi, yRes)
            z = np.linspace(0,np.pi, zRes)

            [X,Y,Z] = np.meshgrid(x,y,z,indexing='ij')

            # Flatten the meshgrid arrays into 1D arrays
            X_flat = X.flatten()
            Y_flat = Y.flatten()
            Z_flat = Z.flatten()

            # Compute the distance from the origin for each point
            distances = np.sqrt(X_flat**2 + Y_flat**2 + Z_flat**2)

            # Filter the points that lie inside the unit sphere
            inside_sphere = distances <= np.pi

            # Create a list of points that are inside the unit sphere
            points_inside_sphere = np.vstack((X_flat[inside_sphere], Y_flat[inside_sphere], Z_flat[inside_sphere])).T

            # Convert to a list of tuples
            points_inside_sphere_list = [tuple(point) for point in points_inside_sphere]

            # print(len(points_inside_sphere_list))
            # fig = plt.figure(figsize=(12, 12))
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(*zip(*points_inside_sphere_list))
            # plt.show()

            orientation_field = np.zeros((4,self.xRes,self.yRes,self.zRes))

            for point in points_inside_sphere_list:
                rotation = sp.spatial.transform.Rotation.from_rotvec(point)
                rotation_matrix = rotation.as_matrix()
                q = (rotation_matrix@self.q.transpose()).transpose()

                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(order_parameter*np.exp(
                            -1j*q[n][0]*self.x - 1j*q[n][1]*self.y)))
                
                Phi2 = np.zeros(self.dims)
                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    Phi2 += abs(eta[n])**2

                theta = np.sqrt(point[0]**2 + point[1]**2 + point[2]**2)
                nx = point[0]/theta if theta != 0 else 0
                ny = point[1]/theta if theta != 0 else 0
                nz = point[2]/theta if theta != 0 else 0

                orientation_field[0] += Phi2*np.cos(2*theta)
                orientation_field[1] += Phi2*np.sin(2*theta)*nx
                orientation_field[2] += Phi2*np.sin(2*theta)*ny
                orientation_field[3] += Phi2*np.sin(2*theta)*nz
            
            return orientation_field
                
    def calc_free_energy(self):
        free_energy_density, _ = self.calc_PFC_free_energy_density_and_chemical_potential()
        return self.calc_integrate_field(free_energy_density)

    #######################################################
    ############### PLOTTING FUNCTIONS ####################
    #######################################################

    def plot_field(self, field, **kwargs):
        """Plots the PFC

        Args:
            field: The field to plot.
            **kwargs: Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.
        Returns:
            ax, fig: The axes and figure containing the plot.
        """
        
        PFC_is_distorted = True if hasattr(self, 'bool_is_distorted') and self.bool_is_distorted else False

        if PFC_is_distorted:
            tool_print_in_color("Note: plotting a strained PFC currently only possible using matplotlib.", 'yellow')
            tool_print_in_color("Output of plot_PFC will therefore be matplotlib ax and fig.", 'yellow')
            kwargs['xmin'] = np.min(self.X)
            kwargs['xmax'] = np.max(self.X)
            kwargs['ymin'] = np.min(self.Y)
            kwargs['ymax'] = np.max(self.Y)
            kwargs['X'] = self.X
            kwargs['Y'] = self.Y
            return plot_field_matplotlib(self, field, **kwargs)
        else:
            return plot_field_plotly(self, field, **kwargs)

    def plot_PFC(self, **kwargs):
        """Plots the PFC

        Args:
            **kwargs: Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.
        Returns:
            ax, fig: The axes and figure containing the plot.
        """
        PFC_has_velocity_field = hasattr(self, 'bool_has_velocity_field') and self.bool_has_velocity_field
    
        if PFC_has_velocity_field:
            return self.plot_field(self.psi[0], **kwargs)
        else:
            return self.plot_field(self.psi, **kwargs)



    def plot_orientation_field(self, orientation_field=None, **kwargs):
        """Plots the orientation field of the phase-field crystal.

        Args:
            orientation_field: The orientation field to plot.
            **kwargs: Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.

        Returns:
            The axes containing the plot.
        """

        if orientation_field is None:
            orientation_field = self.calc_orientation_field()
        
        if self.dim == 2:
            complex_field = orientation_field[0] + 1j*orientation_field[1]
             
            if self.type == 'PhaseFieldCrystal2DTriangular':
                kwargs['cticks'] = [-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]
                kwargs['cticklabels'] = [r'$-\pi/6$', r'$-\pi/9$', r'$-\pi/18$', r'$0$', r'$\pi/18$', r'$\pi/9$', r'$\pi/6$']

            elif self.type == 'PhaseFieldCrystal2DSquare':
                kwargs['cticks'] = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
                kwargs['cticklabels'] = [r'$-\pi/4$', r'$-\pi/8$', r'$0$', r'$\pi/8$', r'$\pi/4$']

            kwargs['plot_method'] = 'phase_angle'

            return self.plot_complex_field(complex_field, **kwargs)

        elif self.dim==3:
            theta = np.arccos(orientation_field[0])/2
            nnorm = np.sqrt(orientation_field[1]**2 + orientation_field[2]**2 + orientation_field[3]**2)
            nx = orientation_field[1]/nnorm
            ny = orientation_field[2]/nnorm 
            nz = orientation_field[3]/nnorm 

            vector_field = theta*np.array([nx,ny,nz])

            return self.plot_vector_field(vector_field, spacing=5, **kwargs)

    def plot_dislocation_nodes(self, dislocation_nodes, **kwargs):
        """
        Plots the dislocation nodes.

        Args:
            dislocation_nodes: The dislocation nodes to plot.
            **kwargs: Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.  
                
        Returns: 
            The axes containing the plot. (matplotlib.axes.Axes)
        """    

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        if self.dim == 2:

            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111)

            x_coords = []
            y_coords = []

            # vx_coords = []
            # vy_coords = []

            Bx_coords = []
            By_coords = []

            for dislocation in dislocation_nodes:

                x_coords.append(dislocation['position'][0])
                y_coords.append(dislocation['position'][1])
                # vx_coords.append(vortex['velocity'][0])
                # vy_coords.append(vortex['velocity'][1])
                Bx_coords.append(dislocation['Burgers_vector'][0])
                By_coords.append(dislocation['Burgers_vector'][1])

            x_coords = np.array(x_coords)
            y_coords = np.array(y_coords)

            # print(x_coords_pos,y_coords_pos)
            # print(x_coords_neg,y_coords_neg)
            ax.scatter(x_coords/self.a0, y_coords/self.a0, marker='o', color='black')

            # ax.quiver(x_coords, y_coords, vx_coords, vy_coords, color='black')
            ax.quiver(x_coords/self.a0, y_coords/self.a0, Bx_coords, By_coords, color='red')
            
            kwargs['ax'] = ax
            self.plot_tool_set_axis_properties(**kwargs)

        elif self.dim == 3:
            # Plotting options
            quiver_scale = 2 # The scale of the quiver arrows

            if ax == None:
                plt.clf()
                ax = plt.gcf().add_subplot(111, projection='3d')

            x_coords = []
            y_coords = []
            z_coords = []

            tx = []
            ty = []
            tz = []

            # vx = []
            # vy = []
            # vz = []

            Bx = []
            By = []
            Bz = []


            for dislocation in dislocation_nodes:
                x_coords.append(dislocation['position'][0])
                y_coords.append(dislocation['position'][1])
                z_coords.append(dislocation['position'][2])

                tx.append(dislocation['tangent_vector'][0])
                ty.append(dislocation['tangent_vector'][1])
                tz.append(dislocation['tangent_vector'][2])

                # vx.append(dislocation['velocity'][0])
                # vy.append(dislocation['velocity'][1])
                # vz.append(dislocation['velocity'][2])

                Bx.append(dislocation['Burgers_vector'][0])
                By.append(dislocation['Burgers_vector'][1])
                Bz.append(dislocation['Burgers_vector'][2])

            tx = np.array(tx)
            ty = np.array(ty)
            tz = np.array(tz)

            # vx = np.array(vx)
            # vy = np.array(vy)
            # vz = np.array(vz)

            Bx = np.array(Bx)
            By = np.array(By)
            Bz = np.array(Bz)

            # if not len(vx) == 0:
            #     v2 =vx**2 + vy**2 + vz**2
            #     v_norm = np.sqrt(max(v2))
            # else:
            #     v_norm = 1

            if not len(Bx) == 0:
                B2 =Bx**2 + By**2 + Bz**2
                B_norm = np.sqrt(max(B2))
            else:
                B_norm = 1

            #ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*tx, quiver_scale*ty, quiver_scale*tz, color='blue')
            # ax.quiver(x_coords, y_coords, z_coords, quiver_scale*vx/v_norm, quiver_scale*vy/v_norm, quiver_scale*vz/v_norm, color='green')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*Bx/B_norm, quiver_scale*By/B_norm, quiver_scale*Bz/B_norm, color='red')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])
            ax.set_zlim([0, self.zmax-self.dz])

            ax.set_aspect('equal')