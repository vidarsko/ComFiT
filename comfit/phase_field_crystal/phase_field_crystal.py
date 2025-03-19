import numpy as np
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint
import copy
from scipy.optimize import curve_fit

# Plot functions
from comfit.plot.plot_field_plotly import plot_field_plotly
from comfit.plot.plot_field_matplotlib import plot_field_matplotlib

from comfit.tool.tool_print_in_color import tool_print_in_color

class PhaseFieldCrystal(BaseSystem):

    def __init__(self, dim, **kwargs):
        """
        Initializes a system to simulate a Phase Field Crystal.
        
        This class is the base of the other phase field crystal models implemented in comfit.
        
        Parameters
        ----------
        dim : int
            The dimension of the system.
        kwargs : dict
            Keyword arguments to set additional parameters. See https://comfitlib.com/ClassPhaseFieldCrystal/
            
        Returns
        -------
        PhaseFieldCrystal
            The system object representing the phase field crystal simulation.
        """
        
        # First initialize the BaseSystem class
        super().__init__(dim, **kwargs)

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
        
        Returns
        -------
        str
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
        """
        Configures the PFC from the amplitudes.
        
        Parameters
        ----------
        eta : array_like, optional
            The amplitudes to configure the PFC from.
        rotation : array_like, optional
            Rotation vector to apply to the crystal.
        
        Returns
        -------
        None
            Configures self.psi and self.psi_f.
        """
        self.psi = self.calc_PFC_from_amplitudes(eta, rotation)
        self.psi_f = self.fft(self.psi)

    def conf_advect_PFC(self, u):
        """
        Advects the PFC according to the displacement field u.
        
        Parameters
        ----------
        u : array_like
            The displacement field to advect the PFC with.
        
        Returns
        -------
        None
            Updates the PFC state.
        """
        self.psi = np.real(self.calc_advect_field(self.psi, u, self.psi_f))
        self.psi_f = self.fft(self.psi)

    def conf_apply_distortion(self, distortion, update_q_and_a_vectors=False):
        """
        Applies a distortion to the PFC.
        
        Parameters
        ----------
        distortion : float or array_like
            The distortion to apply to the PFC. Can be a float for uniform distortion
            or a matrix for more complex distortions.
        update_q_and_a_vectors : bool, optional
            Whether to update the q-vectors and a-vectors, by default False
        
        Returns
        -------
        None
            Updates the PFC state.
        """
        if self.dim == 1:
            self.k[0] = self.k[0]/(1+distortion)
            self.dif[0] = self.dif[0]/(1+distortion)
            self.x = self.x*(1+distortion)
            self.xmax = self.xmax*(1+distortion)
            self.xmin = self.xmin*(1+distortion)
            self.size_x = self.size_x*(1+distortion)
            self.dx = self.dx*(1+distortion)
            self.volume = self.volume*(1+distortion)

        elif self.dim == 2:
            
            distortion_is_shear = False
            if isinstance(distortion, float):
                distortion_matrix = np.array([[1+distortion,0],[0,1+distortion]])
            else:
                distortion_matrix = np.eye(2) + np.array(distortion)

                if abs(distortion_matrix[0,1]) + abs(distortion_matrix[1,0])> 0:
                    distortion_is_shear = True

            # Updating the k and dif vectors
            inverse_distortion_matrix = np.linalg.inv(distortion_matrix)

            if distortion_is_shear:  
                self.bool_is_shear_distorted = True

                # Creating 2D meshgrid
                X,Y = np.meshgrid(self.x.flatten(),self.y.flatten(), indexing='ij')

                # Applying the strain
                X0 = X.copy()
                Y0 = Y.copy()
                X = distortion_matrix[0,0]*X0 + distortion_matrix[0,1]*Y0
                Y = distortion_matrix[1,0]*X0 + distortion_matrix[1,1]*Y0

                # Updating the x and y coordinates
                self.X = X
                self.Y = Y

                k00 = self.k[0].copy()
                k01 = self.k[1].copy()
                self.k[0] = k00*inverse_distortion_matrix[0,0] + k01*inverse_distortion_matrix[0,1]
                self.k[1] = k00*inverse_distortion_matrix[1,0] + k01*inverse_distortion_matrix[1,1]

            else:
                
                x0 = self.x.copy()
                y0 = self.y.copy()
                self.x = x0*distortion_matrix[0,0]
                self.y = y0*distortion_matrix[1,1] 

                X = self.x
                Y = self.y

                k00 = self.k[0].copy()
                k01 = self.k[1].copy()
                self.k[0] = k00*inverse_distortion_matrix[0,0]
                self.k[1] = k01*inverse_distortion_matrix[1,1]

            self.dif[0] = 1j*self.k[0]
            self.dif[1] = 1j*self.k[1]

            if update_q_and_a_vectors:
                    a00 = self.a[:,0].copy()
                    a01 = self.a[:,1].copy()
                    self.a[:,0] = a00*distortion_matrix[0,0] + a01*distortion_matrix[0,1]
                    self.a[:,1] = a00*distortion_matrix[1,0] + a01*distortion_matrix[1,1]

                    q00 = self.q[:,0].copy()
                    q01 = self.q[:,1].copy()
                    self.q[:,0] = q00*inverse_distortion_matrix[0,0] + q01*inverse_distortion_matrix[0,1]
                    self.q[:,1] = q00*inverse_distortion_matrix[1,0] + q01*inverse_distortion_matrix[1,1]

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

            volume_factor = np.linalg.det(distortion_matrix)
            self.dV = self.dV*volume_factor
            self.volume = self.volume*volume_factor


        elif self.dim == 3:
            # Distortion is shear
            distortion_is_shear = False
            if isinstance(distortion, float):
                distortion_matrix = np.array([[1+distortion,0,0],[0,1+distortion,0],[0,0,1+distortion]])
            else:
                distortion_matrix = np.eye(3) + np.array(distortion)

                if abs(distortion_matrix[0,1]) +\
                    abs(distortion_matrix[1,0]) + \
                    abs(distortion_matrix[0,2]) + \
                    abs(distortion_matrix[2,0]) + \
                    abs(distortion_matrix[1,2]) + \
                    abs(distortion_matrix[2,1]) > 0:
                    distortion_is_shear = True

            # Updating the k and dif vectors
            inverse_distortion_matrix = np.linalg.inv(distortion_matrix)

            # if distortion_is_shear:
            if distortion_is_shear:  
                self.bool_is_shear_distorted = True

                # Creating 2D meshgrid
                X,Y,Z = np.meshgrid(self.x.flatten(),self.y.flatten(), self.z.flatten(), indexing='ij')

                # Applying the strain
                X0 = X.copy()
                Y0 = Y.copy()
                Z0 = Z.copy()

                X = distortion_matrix[0,0]*X0 + distortion_matrix[0,1]*Y0 + distortion_matrix[0,2]*Z0
                Y = distortion_matrix[1,0]*X0 + distortion_matrix[1,1]*Y0 + distortion_matrix[1,2]*Z0
                Z = distortion_matrix[2,0]*X0 + distortion_matrix[2,1]*Y0 + distortion_matrix[2,2]*Z0

                # Updating the x and y coordinates
                self.X = X
                self.Y = Y
                self.Z = Z

                k00 = self.k[0].copy()
                k01 = self.k[1].copy()
                k02 = self.k[2].copy()
                self.k[0] = k00*inverse_distortion_matrix[0,0] + k01*inverse_distortion_matrix[0,1] + k02*inverse_distortion_matrix[0,2]
                self.k[1] = k00*inverse_distortion_matrix[1,0] + k01*inverse_distortion_matrix[1,1] + k02*inverse_distortion_matrix[1,2]
                self.k[2] = k00*inverse_distortion_matrix[2,0] + k01*inverse_distortion_matrix[2,1] + k02*inverse_distortion_matrix[2,2]

            else:
                
                x0 = self.x.copy()
                y0 = self.y.copy()
                z0 = self.z.copy()

                self.x = x0*distortion_matrix[0,0]
                self.y = y0*distortion_matrix[1,1] 
                self.z = z0*distortion_matrix[2,2]

                X = self.x
                Y = self.y
                Z = self.z

                k00 = self.k[0].copy()
                k01 = self.k[1].copy()
                k02 = self.k[2].copy()

                self.k[0] = k00*inverse_distortion_matrix[0,0]
                self.k[1] = k01*inverse_distortion_matrix[1,1]
                self.k[2] = k02*inverse_distortion_matrix[2,2]

            self.dif[0] = 1j*self.k[0]
            self.dif[1] = 1j*self.k[1]
            self.dif[2] = 1j*self.k[2]

            if update_q_and_a_vectors:
                    a00 = self.a[:,0].copy()
                    a01 = self.a[:,1].copy()
                    a02 = self.a[:,2].copy()
                    self.a[:,0] = a00*distortion_matrix[0,0] + a01*distortion_matrix[0,1] + a02*distortion_matrix[0,2]
                    self.a[:,1] = a00*distortion_matrix[1,0] + a01*distortion_matrix[1,1] + a02*distortion_matrix[1,2]
                    self.a[:,2] = a00*distortion_matrix[2,0] + a01*distortion_matrix[2,1] + a02*distortion_matrix[2,2]

                    q00 = self.q[:,0].copy()
                    q01 = self.q[:,1].copy()
                    q02 = self.q[:,2].copy()
                    self.q[:,0] = q00*inverse_distortion_matrix[0,0] + q01*inverse_distortion_matrix[0,1] + q02*inverse_distortion_matrix[0,2]
                    self.q[:,1] = q00*inverse_distortion_matrix[1,0] + q01*inverse_distortion_matrix[1,1] + q02*inverse_distortion_matrix[1,2]
                    self.q[:,2] = q00*inverse_distortion_matrix[2,0] + q01*inverse_distortion_matrix[2,1] + q02*inverse_distortion_matrix[2,2]

            # Updating the dx and dy
            original_dx = self.dx
            original_dy = self.dy
            original_dz = self.dz

            self.dx = np.sqrt((distortion_matrix[0,0]*original_dx)**2 + (distortion_matrix[0,1]*original_dy)**2 + (distortion_matrix[0,2]*original_dz)**2)
            self.dy = np.sqrt((distortion_matrix[1,0]*original_dx)**2 + (distortion_matrix[1,1]*original_dy)**2 + (distortion_matrix[1,2]*original_dz)**2)
            self.dz = np.sqrt((distortion_matrix[2,0]*original_dx)**2 + (distortion_matrix[2,1]*original_dy)**2 + (distortion_matrix[2,2]*original_dz)**2)

            # Updating the size_x and size_y
            original_size_x = self.size_x
            original_size_y = self.size_y
            original_size_z = self.size_z

            self.size_x = np.sqrt((distortion_matrix[0,0]*original_size_x)**2 + (distortion_matrix[0,1]*original_size_y)**2 + (distortion_matrix[0,2]*original_size_z)**2)
            self.size_y = np.sqrt((distortion_matrix[1,0]*original_size_x)**2 + (distortion_matrix[1,1]*original_size_y)**2 + (distortion_matrix[1,2]*original_size_z)**2)
            self.size_z = np.sqrt((distortion_matrix[2,0]*original_size_x)**2 + (distortion_matrix[2,1]*original_size_y)**2 + (distortion_matrix[2,2]*original_size_z)**2)

            # Updating the x and y coordinate limits
            self.xmax = np.max(X)
            self.xmin = np.min(X)

            self.ymax = np.max(Y)
            self.ymin = np.min(Y)

            self.zmax = np.max(Z)
            self.zmin = np.min(Z)

            volume_factor = np.linalg.det(distortion_matrix)
            self.dV = self.dV*volume_factor
            self.volume = self.volume*volume_factor

    def conf_strain_to_equilibrium(self):
        """Strain the PFC to equilibrium by adjusting the position variables and k-space variables.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        float
            The final strain value that minimizes the free energy.
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

        pfc_strained = copy.deepcopy(self)
        pfc_strained.conf_apply_distortion(strain, update_q_and_a_vectors=True)

        # Evolve PFC
        pfc_strained.evolve_PFC(number_of_steps, suppress_output=True)

        # Calculate free energy to compare
        average_free_energy_tmp = pfc_strained.calc_free_energy()/pfc_strained.volume

        positive_strain_required = False
        while average_free_energy_tmp < average_free_energy:

            positive_strain_required = True

            average_free_energy = average_free_energy_tmp

            strain += strain_increment

            pfc_strained = copy.deepcopy(self)
            pfc_strained.conf_apply_distortion(strain, update_q_and_a_vectors=True)
            
            pfc_strained.evolve_PFC(number_of_steps, suppress_output=True)

            average_free_energy_tmp = pfc_strained.calc_free_energy()/pfc_strained.volume
            # print('Average free energy at strain: ', strain, ' is: ', average_free_energy_tmp)
        
        if positive_strain_required:
            # Going one back to get the lowest free energy
            final_strain = strain - strain_increment
            self.conf_apply_distortion(final_strain, update_q_and_a_vectors=True)
            # tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')

        else: #negative strain required

            strain = - strain_increment

            pfc_strained = copy.deepcopy(self)
            pfc_strained.conf_apply_distortion(strain, update_q_and_a_vectors=True)

            pfc_strained.evolve_PFC(number_of_steps, suppress_output=True)
            average_free_energy_tmp = pfc_strained.calc_free_energy()/pfc_strained.volume
            # print('Average Free energy at strain: ', strain, ' is: ', average_free_energy_tmp)

            while average_free_energy_tmp < average_free_energy:

                average_free_energy = average_free_energy_tmp

                strain -= strain_increment

                pfc_strained = copy.deepcopy(self)
                pfc_strained.conf_apply_distortion(strain, update_q_and_a_vectors=True)

                pfc_strained.evolve_PFC(number_of_steps, suppress_output=True)

                average_free_energy_tmp = pfc_strained.calc_free_energy()/pfc_strained.volume

                # print('Average free energy at strain: ', strain, ' is: ', average_free_energy_tmp)

            # Going one back to get the lowest free energy
            final_strain = strain + strain_increment
            self.conf_apply_distortion(final_strain, update_q_and_a_vectors=True)
            # tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')   

        return final_strain

    def conf_create_polycrystal(self, type, **kwargs):
        """Creates a polycrystal.
        
        Parameters
        ----------
        type : str
            The type of polycrystal to create ('circular' or 'four_grain').
        kwargs : dict
            Additional arguments for the polycrystal creation, including:
            
            - relaxation_time : float
                The relaxation time to use for the polycrystal creation.
            - rotation : float
                The rotation angle for 'circular' type (default: pi/6).
            - position : array_like
                The position for 'circular' type (default: system midpoint).
            - radius : float
                The radius for 'circular' type (default: size_min/4).
        
        Returns
        -------
        None
            Updates the PFC state.
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
            self.psi_f = self.fft(self.psi)

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
            pfcRotated = np.roll(pfcRotated, -round(self.yRes/2), axis=1)
            pfcRotated = np.roll(pfcRotated, -round(self.xRes/3), axis=0)
            region = np.bool_((l4*~(l7)*~(l8) + ~(l1)*~(l3)*~(l7))*zdir)
            self.psi[region] = pfcRotated[region]

            pfcRotated = self.calc_PFC_from_amplitudes(self.eta0, rotation=[0,0,22.5/180*np.pi])
            pfcRotated = np.roll(pfcRotated, -round(self.yRes/2), axis=1)
            pfcRotated = np.roll(pfcRotated, round(self.xRes/3), axis=0)
            region = np.bool_((l1*l6*l7 + l7*l10*~(l4))*zdir)
            self.psi[region] = pfcRotated[region]

            pfcRotated = self.calc_PFC_from_amplitudes(self.eta0, rotation=[0,0,45/180*np.pi])
            pfcRotated = np.roll(pfcRotated, -round(self.xRes/2), axis=0)
            region = np.bool_((l1*~(l4)*~(l5) + ~(l1)*l4*l9)*zdir)
            self.psi[region] = pfcRotated[region]

            self.psi_f = self.fft(self.psi)

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
        
        Parameters
        ----------
        number_of_steps : int
            The number of steps to evolve the PFC.
        method : str, optional
            The method to use for the evolution, by default 'ETD2RK'.
        suppress_output : bool, optional
            Whether to suppress the progress bar, by default False.
        
        Returns
        -------
        None
            Updates self.psi and self.psi_f.
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
                                            non_linear_evolution_function_f,
                                            self.psi, self.psi_f)

            # These steps seem to be necessary for numerical stability (Vidar 18.12.23)
            self.psi = np.real(self.psi)
            self.psi_f = self.fft(self.psi)

    # Nonlinear part conserved
    def calc_nonlinear_evolution_function_conserved_f(self, psi, t):
        """Calculate the nonlinear part of the evolution function for conserved dynamics.
        
        Parameters
        ----------
        psi : ndarray
            The phase field.
        t : float
            The time.
            
        Returns
        -------
        ndarray
            The nonlinear part of the evolution function in Fourier space.
        """
        return -self.calc_k2()*self.fft(self.t * psi ** 2 + self.v * psi ** 3)

    # Non-linear part unconserved
    def calc_nonlinear_evolution_function_unconserved_f(self, psi, t):
        """Calculate the nonlinear part of the evolution function for unconserved dynamics.
        
        Parameters
        ----------
        psi : ndarray
            The phase field.
        t : float
            The time.
            
        Returns
        -------
        ndarray
            The nonlinear part of the evolution function in Fourier space.
        """
        return -self.fft(self.t * psi ** 2 + self.v * psi ** 3)
    #######################################################

    #######################################################
    ##### MECHANICAL EQUILIBRIUM (CONSERVED) ##############
    #######################################################
    def evolve_PFC_mechanical_equilibrium(self, time, Delta_t = 10, method='ETD2RK'):
        """Evolves the PFC in mechanical equilibrium.
        
        Parameters
        ----------
        time : float
            The total time to evolve the PFC.
        Delta_t : float, optional
            The time step for the mechanical equilibrium evolution, by default 10.
        method : str, optional
            The method to use for the evolution, by default 'ETD2RK'.
            
        Returns
        -------
        None
            Updates self.psi and self.psi_f.
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
                                gamma_S = 2**-6,
                                rho0 = 2**-6):
        """Evolves the PFC according to hydrodynamic PFC dynamics.
        
        This requires introducing a velocity field. 
        If psi does not contain this field, it is added to the components psi[1], psi[2], psi[3].
        
        Parameters
        ----------
        number_of_steps : int
            The number of steps to evolve the PFC.
        method : str, optional
            The method to use for the evolution, by default 'ETD2RK'.
        gamma_S : float, optional
            The surface tension coefficient, by default 2**-6.
        rho0 : float, optional
            The mass density, by default 2**-6.
            
        Returns
        -------
        None
            Updates self.psi and self.psi_f.
        """
        
        if hasattr(self,'bool_has_velocity_field'):
            pass
        else:
            self.bool_has_velocity_field = True
            self.psi = np.array([self.psi]+[np.zeros_like(self.psi)]*self.dim)
            self.psi_f = np.array([self.psi_f]+[np.zeros_like(self.psi_f)]*self.dim)
            # print("psi shape", self.psi.shape).

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
            self.psi_f = self.fft(self.psi )

    # Linear part
    def calc_omega_hydrodynamic_f(self):
        """Calculates the hydrodynamic evolution function omega_f.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        ndarray
            The hydrodynamic evolution function omega_f in Fourier space.
        """
        k2 = self.calc_k2()

        if self.type_of_evolution == 'conserved':
            omega_f = -self.calc_k2()*(self.r + self.calc_L_f()**2)
        elif self.type_of_evolution == 'unconserved':
            omega_f = -(self.r + self.calc_L_f()**2)

        return np.array([omega_f]+[-self.gamma_S/self.rho0*k2]*self.dim)

    # Nonlinear part
    def calc_nonlinear_hydrodynamic_evolution_function_f(self, field, t):
        """Calculates the hydrodynamic evolution function of the PFC.
        
        Parameters
        ----------
        field : ndarray
            The field to calculate the evolution function of.
        t : float
            The time.
            
        Returns
        -------
        ndarray
            The nonlinear evolution function for the hydrodynamic PFC in Fourier space.
        """

        field_f = self.fft(field)

        k2 = self.calc_k2()

        N0_f = -k2*self.fft(self.t * field[0] ** 2 + self.v * field[0] ** 3) \
            - self.fft(sum([field[i+1]*self.ifft(self.dif[i]*field_f[0]) for i in range(self.dim)]))
        
        force_density_f = self.calc_stress_divergence_f(field_f[0])

        return np.array([N0_f] + [1/self.rho0*(force_density_f[i]+self.external_force_density_f[i]) for i in range(self.dim)])
    #######################################################


    #######################################################
    #######################################################
    ############## CALCULATION FUNCTIONS ##################
    #######################################################
    #######################################################
    
    def calc_strained_amplitudes(self):
        """Strains the PFC to equilibrium and returns the amplitudes.
    
        Parameters
        ----------
        None

        Returns
        -------
        tuple
            Depending on the PFC type, returns either:
            - (final_strain, psi0, A, el_lambda, el_mu, el_gamma) for 1 independent amplitude
            - (final_strain, psi0, A, B, el_lambda, el_mu, el_gamma) for 2 independent amplitudes
            - (final_strain, psi0, A, B, C, el_lambda, el_mu, el_gamma) for 3 independent amplitudes
        """
        # tool_print_in_color('Proto amplitudes and elastic constants', 'blue')
        # tool_print_in_color('---', 'blue')
        # print(f'Proto psi0: {self.psi0:.02f}')
        # # print(f'Proto A: {self.A:.05f}')
        # if self.type in ['PhaseFieldCrystal2DSquare','PhaseFieldCrystal3DFaceCenteredCubic','PhaseFieldCrystal3DSimpleCubic']:
        #     print(f'Proto B: {self.B:.05f}')
        # if self.type in ['PhaseFieldCrystal3DSimpleCubic']:
        #     print(f'Proto C: {self.C:.05f}')

        # print(f'Proto mu: {self.el_mu:.05f}')
        # print(f'Proto lambda: {self.el_lambda:.05f}')
        # print(f'Proto gamma: {self.el_gamma:.05f}')

        print('---')
        tool_print_in_color('Straining PFC to reach equilibrium', 'blue')
        print('---')

        # Strain PFC to equilibrium
        self.conf_PFC_from_amplitudes()
        final_strain = self.conf_strain_to_equilibrium()
        self.evolve_PFC(200, suppress_output=True)

        # tool_print_in_color('Amplitudes after strain', 'blue')
        # tool_print_in_color('---', 'blue')
        print(f'Equilibrium strain: {final_strain:.05f}')
        print(f'Equilibrium q-vector: {1/(1+final_strain):.05f}')

        psi0 = np.mean(self.psi)
        # print(f'Eq. psi0: {psi0:.02f}')

        eta = self.calc_demodulate_PFC(only_primary_modes=False)
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

        print('---')
        tool_print_in_color('Equilibrium amplitudes:', 'blue')
        print('---')
        print(f'Proto mean density psi0      : {self.psi0:.02f}')
        print(f'Equilibrium mean density psi0: {psi0:.02f}')
        print('Ratio (equilibrium/proto)     : {:.05f}'.format(psi0/self.psi0))

        if number_of_independent_amplitudes >= 1:
            print('---')
            print(f'Proto amplitude       A : {self.A:.05f}')
            print(f'Equilibrium amplitude A : {A:.05f}')
            print('Ratio (equilibrium/proto): {:.05f}'.format(A/self.A))
        if number_of_independent_amplitudes >= 2:
            print('---')
            print(f'Proto amplitude       B : {self.B:.05f}')
            print(f'Equilibrium amplitude B : {B:.05f}')
            print('Ratio (equilibrium/proto): {:.05f}'.format(B/self.B))
        
        if number_of_independent_amplitudes >= 3:
            print('---')
            print(f'Proto amplitude       C : {self.C:.05f}')
            print(f'Equilibrium amplitude C : {C:.05f}')
            print('Ratio (equilibrium/proto): {:.05f}'.format(C/self.C))
        
        # Finding elastic constants
        # def elastic_energy(strain, el_lambda, el_mu, el_gamma):
        #     exx, exy, eyy = strain
        #     return el_lambda/2*(exx+eyy)**2 + el_mu*(exx**2 + 2*exy**2 + eyy**2) + el_gamma/2*(exx**2 + eyy**2)

        # strain_magnitudes = np.linspace(-0.01,0.01,21)
        # exx = np.array([[a,0,a] for a in strain_magnitudes]).flatten()
        # exy = np.array([[0,a,0] for a in strain_magnitudes]).flatten()
        # eyy = np.array([[a,0,-a] for a in strain_magnitudes]).flatten()

        f0 = self.calc_free_energy()/self.volume

        # free_energies = np.zeros_like(exx)
        # for n in range(len(exx)):
        #     if self.dim == 2:
        #         distortion = [[exx[n], exy[n]],[exy[n], eyy[n]]]
        #     elif self.dim == 3:
        #         distortion = [[exx[n], exy[n], 0],[exy[n], eyy[n], 0],[0,0,0]]

        #     pfc_strained = copy.deepcopy(self)
        #     pfc_strained.conf_apply_distortion(distortion)
        #     f = pfc_strained.calc_free_energy()/pfc_strained.volume
        #     free_energies[n] = f-f0

        # params,_ = curve_fit(elastic_energy, (exx, exy, eyy), free_energies)
        # el_lambda, el_mu, el_gamma = params

        epsilon=0.01
        
        if self.dim == 2:
            pure_shear_distortion = [[0,epsilon],[epsilon,0]]
            volume_conserving_compression = [[epsilon,0],[0,-epsilon]]
            pure_compression_pos = [[epsilon,0],[0,epsilon]]
            pure_compression_neg = [[-epsilon,0],[0,-epsilon]]
        else:
            pure_shear_distortion = [[0,epsilon,0],[epsilon,0,0],[0,0,0]]
            volume_conserving_compression = [[epsilon,0,0],[0,-epsilon,0],[0,0,0]]
            pure_compression_pos = [[epsilon,0,0],[0,epsilon,0],[0,0,0]]
            pure_compression_neg = [[-epsilon,0,0],[0,-epsilon,0],[0,0,0]]

        pfc_strained = copy.deepcopy(self)
        pfc_strained.conf_apply_distortion(pure_shear_distortion)
        f_pure_shear = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        pfc_strained = copy.deepcopy(self)
        pfc_strained.conf_apply_distortion(volume_conserving_compression)
        f_volume_conserving_compression = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        pfc_strained = copy.deepcopy(self)
        pfc_strained.conf_apply_distortion(pure_compression_pos)
        f_pure_compression_pos = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        pfc_strained = copy.deepcopy(self)
        pfc_strained.conf_apply_distortion(pure_compression_neg)
        f_pure_compression_neg = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

        f_pure_compression = (f_pure_compression_pos + f_pure_compression_neg)/2

        el_lambda = (0.5*f_pure_compression - 0.5*f_volume_conserving_compression)/(epsilon**2)
        el_mu = (0.5*f_pure_shear)/(epsilon**2)
        el_gamma = (-f_pure_shear + f_volume_conserving_compression)/(epsilon**2)


        if self.type == 'PhaseFieldCrystal1DPeriodic':
            el_lambda_from_eq_amplitudes = 2*(A**2)
        elif self.type == 'PhaseFieldCrystal2DTriangular':
            el_lambda_from_eq_amplitudes = 3*(A**2)
            el_mu_from_eq_amplitudes = 3*(A**2)
            el_gamma_from_eq_amplitudes = 0
        elif self.type == 'PhaseFieldCrystal2DSquare':
            el_lambda_from_eq_amplitudes = 16*(B**2)
            el_mu_from_eq_amplitudes = 16*(B**2)
            el_gamma_from_eq_amplitudes = 8*A**2 - 32*B**2
        elif self.type == 'PhaseFieldCrystal3DBodyCenteredCubic':
            el_lambda_from_eq_amplitudes = 4*(A**2)
            el_mu_from_eq_amplitudes = 4*(A**2)
            el_gamma_from_eq_amplitudes = -4*(A**2)
        elif self.type == 'PhaseFieldCrystal3DFaceCenteredCubic':
            el_lambda_from_eq_amplitudes = 32/81*(A**2)
            el_mu_from_eq_amplitudes = 32/81*(A**2)
            el_gamma_from_eq_amplitudes = 32/81*(2*B**2-A**2)
        elif self.type == 'PhaseFieldCrystal3DSimpleCubic':
            el_lambda_from_eq_amplitudes = 16*(B**2) + 128*(C**2)
            el_mu_from_eq_amplitudes = 16*(B**2) + 128*(C**2)
            el_gamma_from_eq_amplitudes = 32*(A**2) - 16*(B**2) - 256*(C**2)

        print('---')
        tool_print_in_color('Equilibrium elastic constants:', 'blue')
        print('---')
        print('Proto lambda                             : {:.05f}'.format(self.el_lambda))
        print('Equilibrium lambda (numerical)           : {:.05f}'.format(el_lambda))
        print('Equilibrium lambda (from eq. amplitudes) : {:.05f}'.format(el_lambda_from_eq_amplitudes))
        print('Ratio (equilibrium/proto)                : {:.05f}'.format(el_lambda/self.el_lambda))

        print('---')
        print('Proto mu                             : {:.05f}'.format(self.el_mu))
        print('Equilibrium mu (numerical)           : {:.05f}'.format(el_mu))
        print('Equilibrium mu (from eq. amplitudes) : {:.05f}'.format(el_mu_from_eq_amplitudes))
        print('Ratio (equilibrium/proto)            : {:.05f}'.format(el_mu/self.el_mu))

        print('---')
        print('Proto gamma                             : {:.05f}'.format(self.el_gamma))
        print('Equilibrium gamma (numerical)           : {:.05f}'.format(el_gamma))
        print('Equilibrium gamma (from eq. amplitudes) : {:.05f}'.format(el_gamma_from_eq_amplitudes))
        print('Ratio (equilibrium/proto)               : {:.05f}'.format(el_gamma/self.el_gamma))

        if number_of_independent_amplitudes == 1:
            return final_strain, psi0, A, el_lambda, el_mu, el_gamma
        elif number_of_independent_amplitudes == 2:
            return final_strain, psi0, A, B, el_lambda, el_mu, el_gamma
        elif number_of_independent_amplitudes == 3:
            return final_strain, psi0, A, B, C, el_lambda, el_mu, el_gamma

    def calc_PFC_free_energy_density_and_chemical_potential(self, field=None, field_f=None):
        """Calculates the free energy density and chemical potential of the PFC.

        Parameters
        ----------
        field : ndarray, optional
            The field to calculate the free energy density and chemical potential of.
            If None, self.psi is used.
        field_f : ndarray, optional
            The Fourier transform of the field.
            If None, self.psi_f is used.

        Returns
        -------
        tuple
            A tuple containing:
            - free_energy_density (ndarray): The free energy density of the PFC
            - chemical_potential (ndarray): The chemical potential of the PFC
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

        Lpsi = self.ifft(self.calc_L_f()*psi_f)

        free_energy_density = 1/2*Lpsi**2 \
            + 1/2*self.r*psi2 + 1/3*self.t*psi3 + 1/4*self.v*psi4
        
        # print("Lpsi shape",self.calc_Lpsi(psi_f).shape)

        L2psi = self.ifft(self.calc_L_f()**2*psi_f)

        chemical_potential = L2psi*psi \
            + self.r*psi + self.t*psi2 + self.v*psi3
        
        # print("L2psi shape",self.calc_L2psi(psi_f).shape)
        
        # print("Free energy density shape", free_energy_density.shape)
        # print("Chem pot shape", chemical_potential.shape)
        return np.real(free_energy_density), np.real(chemical_potential)
    

    def calc_displacement_field_to_equilibrium(self):
        """Calculates the displacement field needed to put the PFC in mechanical equilibrium.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            The displacement field u that would bring the system to mechanical equilibrium
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

        return np.real(self.ifft(u_f))

    # Initial configuration methods
    def calc_amplitudes_with_dislocation(self, eta=None, x=None, y=None, dislocation_type=1):
        """Calculate the amplitudes with a single point dislocation inserted.

        Parameters
        ----------
        eta : array_like, optional
            The amplitudes to insert the dislocation in. If None, uses default amplitudes.
        x : float, optional
            The x-coordinate of the dislocation. If None, uses system midpoint.
        y : float, optional
            The y-coordinate of the dislocation. If None, uses system midpoint.
        dislocation_type : int, optional
            The dislocation type to insert, default is 1.
        
        Returns
        -------
        ndarray
            The amplitudes containing the dislocation.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 2.
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
        """Insert a dislocation dipole in the system corresponding to dislocation type and its negative.

        Parameters
        ----------
        eta : array_like, optional
            The amplitudes to insert the dislocation dipole in. If None, uses default amplitudes.
        x1 : float, optional
            The x-coordinate of the first dislocation. If None, uses 1/3 of system width.
        y1 : float, optional
            The y-coordinate of the first dislocation. If None, uses 1/2 of system height.
        x2 : float, optional
            The x-coordinate of the second dislocation. If None, uses 2/3 of system width.
        y2 : float, optional
            The y-coordinate of the second dislocation. If None, uses 1/2 of system height.
        dislocation_type : int, optional
            The dislocation type to insert, default is 1.
        
        Returns
        -------
        ndarray
            The amplitudes with the dislocation dipole inserted.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 2.
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

        if eta is None:
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
        """Insert a dislocation ring in the system corresponding to dislocation type.

        Parameters
        ----------
        eta : array_like, optional
            The amplitudes to insert the dislocation ring in. If None, uses default amplitudes.
        position : array_like, optional
            The position of the dislocation ring. If None, uses system midpoint.
        radius : float, optional
            The radius of the dislocation ring. If None, uses 1/3 of minimum system dimension.
        normal_vector : array_like, optional
            The normal vector of the dislocation ring, default is [0, 0, 1].
        dislocation_type : int, optional
            The dislocation type to insert, default is 1.

        Returns
        -------
        ndarray
            The amplitudes with the dislocation ring inserted.
            
        Raises
        ------
        Exception
            If the dimension of the system is not 3.
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
        """Calculate the PFC from the amplitudes.

        Parameters
        ----------
        eta : array_like, optional
            The amplitudes to calculate the PFC from. If None, uses default amplitudes.
        rotation : array_like or float, optional
            The rotation to apply to the crystal. If float, rotation around z-axis.
            If array, rotation vector [rx, ry, rz].

        Returns
        -------
        ndarray
            The phase field crystal density field.
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

    def calc_demodulate_PFC(self, only_primary_modes=True):
        """Demodulate the PFC to extract amplitudes.

        Parameters
        ----------
        only_primary_modes : bool, optional
            Whether to only extract primary reciprocal lattice modes, default is True.

        Returns
        -------
        ndarray
            The amplitudes corresponding to the demodulated PFC.
        """
        
        number_of_modes = self.number_of_primary_reciprocal_lattice_modes if only_primary_modes else self.number_of_reciprocal_lattice_modes

        eta = np.zeros([number_of_modes] + self.dims, 
                       dtype=complex)

        Gaussian_filter_f = self.calc_Gaussian_filter_f()

        order_parameter = self.psi if self.psi.ndim == self.dim else self.psi[0]

        if self.dim == 2:
                for n in range(number_of_modes):
                    eta[n] = self.ifft(Gaussian_filter_f*self.fft(order_parameter*np.exp(
                        -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y)))

        elif self.dim == 3:
            for n in range(number_of_modes):
                eta[n] = self.ifft(Gaussian_filter_f*self.fft(order_parameter*np.exp(
                    -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y - 1j*self.q[n][2]*self.z  
                    )))
                
        return eta

    def calc_stress_tensor_microscopic(self):
        """Calculate the microscopic stress of the phase-field crystal.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            The microscopic stress tensor of the phase-field crystal.
            For 2D: tensor with components [σxx, σxy, σyy].
            For 3D: tensor with components [σxx, σxy, σxz, σyy, σyz, σzz].
            
        Raises
        ------
        Exception
            If the dimension of the system is 1D.
        """
        if self.dim==1:
            raise Exception("The stress tensor is not yet defined in 1D.")
        elif self.dim==2:
            stress = np.zeros((3,self.xRes,self.yRes))

            Lpsi = np.real(self.ifft(self.calc_L_f()*self.psi_f))
            stress[0] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[0]*self.dif[0]*self.psi_f))
            stress[1] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[0]*self.dif[1]*self.psi_f))
            stress[2] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[1]*self.dif[1]*self.psi_f))
        
        elif self.dim==3:
            stress = np.zeros((6,self.xRes,self.yRes,self.zRes))

            Lpsi = np.real(self.ifft(self.calc_L_f()*self.psi_f))
            stress[0] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[0]*self.dif[0]*self.psi_f))
            stress[1] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[0]*self.dif[1]*self.psi_f))
            stress[2] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[0]*self.dif[2]*self.psi_f))
            stress[3] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[1]*self.dif[1]*self.psi_f))
            stress[4] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[1]*self.dif[2]*self.psi_f))
            stress[5] = -2*Lpsi*np.real(self.ifft(self.calc_L_sum_f()*self.dif[2]*self.dif[2]*self.psi_f))

        return stress


    def calc_stress_tensor(self):
        """Calculates the stress of the phase-field crystal.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            The stress tensor.
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
            stress[n] = Gaussian_filter_f*self.fft(stress[n])

        # return stress
        # return np.real(self.ifft(Gaussian_filter_f))
        return np.real(self.ifft(stress))

    def calc_stress_divergence_f(self, field_f = None):
        """Calculates the divergence of the stress tensor in Fourier space.

        Parameters
        ----------
        field_f : ndarray, optional
            The field in Fourier space. If None, uses the current system field.

        Returns
        -------
        ndarray
            The divergence of the stress tensor in Fourier space.
        """
        if field_f is None:
            PFC_has_velocity_field = hasattr(self, 'bool_has_velocity_field') and self.bool_has_velocity_field
            if PFC_has_velocity_field:
                field_f = self.psi_f[0]
            else:
                field_f = self.psi_f

        L_f = self.calc_L_f()
        Lpsi = self.ifft(L_f*field_f)

        L_sum_f = self.calc_L_sum_f()
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*self.fft(
                sum([
                self.ifft(L_f*self.dif[i]*field_f)*self.ifft(L_sum_f*self.dif[i]*self.dif[j]*field_f) 
                for i in range(self.dim)
                ]) 
                +Lpsi*self.ifft(L_sum_f*self.dif[j]*(-k2)*field_f)) 
                for j in range(self.dim)]
                )


    def calc_structure_tensor_f(self):
        """Calculates the structure tensor of the phase-field crystal in Fourier space.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            The structure tensor in Fourier space.
        """

        if hasattr(self,'bool_has_velocity_field') and self.bool_has_velocity_field:
            field_f = self.psi_f[0]
        else:
            field_f = self.psi_f

        # Calculate the gradient
        diPsi = np.zeros([self.dim] + self.dims, dtype=complex)
        for i in range(self.dim):
            diPsi[i] = self.dif[i]*field_f
        diPsi = np.real(self.ifft(diPsi))

        if self.dim == 1:
            number_of_independent_strain_components = 1
        elif self.dim == 2:
            number_of_independent_strain_components = 3
        elif self.dim == 3:
            number_of_independent_strain_components = 6

        structure_tensor_f = np.zeros([number_of_independent_strain_components] + self.dims, dtype=complex)

        Gaussian_filter_f = self.calc_Gaussian_filter_f()
        if self.dim == 1:
            structure_tensor_f[0] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[0])
        elif self.dim == 2:
            structure_tensor_f[0] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[0])
            structure_tensor_f[1] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[1])
            structure_tensor_f[2] = Gaussian_filter_f*self.fft(diPsi[1]*diPsi[1])
        elif self.dim == 3:
            structure_tensor_f[0] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[0])
            structure_tensor_f[1] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[1])
            structure_tensor_f[2] = Gaussian_filter_f*self.fft(diPsi[0]*diPsi[2])
            structure_tensor_f[3] = Gaussian_filter_f*self.fft(diPsi[1]*diPsi[1])
            structure_tensor_f[4] = Gaussian_filter_f*self.fft(diPsi[1]*diPsi[2])
            structure_tensor_f[5] = Gaussian_filter_f*self.fft(diPsi[2]*diPsi[2])

        return structure_tensor_f
            

    def calc_structure_tensor(self):
        """Calculates the structure tensor of the phase-field crystal.

        Parameters
        ----------
        None
        
        Returns
        -------
        ndarray
            The structure tensor.
        """
        structure_tensor_f = self.calc_structure_tensor_f()
        return np.real(self.ifft(structure_tensor_f))
        

    def calc_strain_tensor(self):
        """Calculates the strain of the phase-field crystal.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            The strain tensor.
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
        """Calculates the dislocation density.

        Parameters
        ----------
        eta : ndarray, optional
            The amplitudes to calculate the dislocation density from.
            If None, amplitudes are calculated using demodulation.

        Returns
        -------
        ndarray
            The dislocation density.
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
        """Calculates the dislocation nodes.

        Parameters
        ----------
        None
        
        Returns
        -------
        list
            A list of dictionaries containing information about each dislocation node.
            
        Raises
        ------
        Exception
            If the PFC is distorted.
        """

        # At the moment (24.09.24 - Vidar), this function only works for unstrained PFCs
        if hasattr(self,'bool_is_shear_distorted') and self.bool_is_shear_distorted:
            raise Exception("Dislocation nodes cannot be calculated for distorted PFCs.")

        alpha = self.calc_dislocation_density()

        if self.dim == 2:
            # self.plot_field(np.sqrt(alpha[0]**2 + alpha[1]**2))
            # plt.show()
            # print(alpha)
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
                    eta[n] = self.ifft(Gaussian_filter_f*self.fft(order_parameter*np.exp(
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
                    eta[n] = self.ifft(Gaussian_filter_f*self.fft(order_parameter*np.exp(
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
        """
        Calculate the total free energy of the system by integrating the free energy density over the computational domain.
        
        This method computes the phase field crystal free energy density using `calc_PFC_free_energy_density_and_chemical_potential`,
        and then integrates it over the entire domain using `calc_integrate_field`.
        
        Returns
        -------
        float
            The total free energy of the system.
        
        See Also
        --------
        calc_PFC_free_energy_density_and_chemical_potential : Method that calculates the free energy density and chemical potential.
        calc_integrate_field : Method that integrates a field over the computational domain.
        """
        free_energy_density, _ = self.calc_PFC_free_energy_density_and_chemical_potential()
        return self.calc_integrate_field(free_energy_density)

    #######################################################
    ############### PLOTTING FUNCTIONS ####################
    #######################################################
    def plot_field(self, field, **kwargs):
        """Plots the PFC.

        Parameters
        ----------
        field : ndarray
            The field to plot.
        kwargs : dict
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.
            
        Returns
        -------
        tuple
            A tuple containing (ax, fig), the axes and figure containing the plot.
        """
        
        PFC_is_distorted = True if hasattr(self, 'bool_is_shear_distorted') and self.bool_is_shear_distorted else False

        if PFC_is_distorted:
            kwargs['xmin'] = np.min(self.X)
            kwargs['xmax'] = np.max(self.X)
            kwargs['ymin'] = np.min(self.Y)
            kwargs['ymax'] = np.max(self.Y)
            kwargs['X'] = self.X
            kwargs['Y'] = self.Y
        
        if self.plot_lib == 'matplotlib':
            return plot_field_matplotlib(self, field, **kwargs)
        elif self.plot_lib == 'plotly':
            return plot_field_plotly(self, field, **kwargs)

    def plot_PFC(self, **kwargs):
        """Plots the PFC.

        Parameters
        ----------
        kwargs : dict
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.
            
        Returns
        -------
        tuple
            A tuple containing (ax, fig), the axes and figure containing the plot.
        """
        PFC_has_velocity_field = hasattr(self, 'bool_has_velocity_field') and self.bool_has_velocity_field

        kwargs['colormap'] = kwargs.get('colormap', 'viridis' if self.plot_lib == 'plotly' else 'viridis')
        if PFC_has_velocity_field:
            return self.plot_field(self.psi[0], **kwargs)
        else:
            return self.plot_field(self.psi, **kwargs)

    def plot_orientation_field(self, orientation_field=None, **kwargs):
        """Plots the orientation field of the phase-field crystal.

        Parameters
        ----------
        orientation_field : ndarray, optional
            The orientation field to plot. If None, it will be calculated.
        kwargs : dict
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.

        Returns
        -------
        tuple
            A tuple containing (ax, fig), the axes and figure containing the plot.
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
