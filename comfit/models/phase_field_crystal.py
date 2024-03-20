import numpy as np
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint


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

        self.Phi = 2*sum(np.array(self.eta0)**2)


    def __str__(self):
        return self.type

    # CONFIGURATION FUNCTIONS
    def conf_PFC_from_amplitudes(self, eta=None):

        self.psi = self.psi0

        if eta is None:
            eta = self.eta0

        for n in range(self.number_of_reciprocal_lattice_modes):

            if self.dim == 1:
                self.psi += 2 * eta[n] * np.exp(1j * self.q[n][0] * self.x)

            elif self.dim == 2:
                self.psi += 2 * eta[n] * np.exp(1j * (self.q[n][0] * self.x + self.q[n][1] * self.y))

            elif self.dim == 3:
                self.psi += 2 * eta[n] * np.exp(1j * (self.q[n][0] * self.x + self.q[n][1] * self.y + self.q[n][2] * self.z))

        self.psi = np.real(self.psi)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_advect_PFC(self,u):
        """
        Advects the PFC according to the displacement field u.
        """

        self.psi = np.real(self.calc_advect_field(self.psi, u, self.psi_f))
        self.psi_f = sp.fft.fftn(self.psi)

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
        self.psi_f = sp.fft.fftn(self.psi)

    def evolve_PFC_mechanical_equilibrium(self, time, Delta_t = 10, method='ETD2RK'):
        """
        Evolves the PFC in mechanical equilibrium. 

        Input:
            number_of_steps: The number of steps to evolve the PFC
            Delta_t: The time step for the mechanical equilibrium evolution
            method: The method to use for the evolution (default: ETD2RK)
        """
        
        number_of_steps = round(time/self.dt)
        number_of_steps_per_iteration = round(Delta_t/self.dt)
        number_of_iterations = round(number_of_steps/number_of_steps_per_iteration)

        for n in tqdm(range(number_of_iterations), desc='Evolving the PFC in mechanical equilibrium'):
            self.conf_advect_PFC(self.calc_displacement_field_to_equilibrium())
            self.evolve_PFC(number_of_steps_per_iteration, method)
            


    def evolve_PFC_hydrodynamic(self, number_of_steps, 
                                method = 'ETD2RK',
                                gamma_S = 2**-2,
                                rho0 = 2**-2):
        """
        We extend psi to also contain the hydrodynamic field in components psi[1],psi[2],psi[3]
        """
        
        if hasattr(self,'velocity_field'):
            pass
        else:
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
        self.psi_f = sp.fft.fftn(self.psi, axes = (range ( - self.dim , 0) ) )


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

    def calc_nonlinear_hydrodynamic_evolution_function_f(self, field, t):

        field_f = sp.fft.fftn(field, axes =( range ( - self . dim , 0) ))

        k2 = self.calc_k2()

        N0_f = -k2*sp.fft.fftn(self.t * field[0] ** 2 + self.v * field[0] ** 3) \
            - sp.fft.fftn(sum([field[i+1]*sp.fft.ifftn(self.dif[i]*field_f[0]) for i in range(self.dim)]))
        
        force_density_f = self.calc_stress_divergence_f(field_f[0])

        return np.array([N0_f] + [1/self.rho0*force_density_f[i] for i in range(self.dim)])

    def calc_nonlinear_evolution_function_f(self, psi, t):
        return -self.calc_k2()*sp.fft.fftn(self.t * psi ** 2 + self.v * psi ** 3)

    def calc_displacement_field_to_equilibrium(self):
        """
        Calculates the displacement field needed to put the PFC in mechanical equilibrium.

        Input:
            None
        
        Output:
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

        :param eta: The amplitudes currently being manipulated (default: self.eta0)
        :param x1: The x-position of the first dislocation (default: 1/3 xmax)
        :param y1: The y-position of the first dislocation (default: 1/2 ymax)
        :param x2: The x-position of the second dislocation (default: 2/3 xmax)
        :param y2: The y-position of the second dislocation (default: 1/2 ymax)
        :param dislocation_type: the type of dislocation dipole (default: 1)
        :Output:
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
        """
        Calculates the PFC from the amplitudes.
        """ 

        psi = self.psi0

        if rotation is None:
            rotation = [0,0,0]
        else:
            if self.dim == 2:
                rotation = [0,0,rotation]

        rotation = sp.spatial.transform.Rotation.from_rotvec(rotation)
        rotation_matrix = rotation.as_matrix()

        # Rotate q-vectors to the new orientation
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
        eta = np.zeros([self.number_of_primary_reciprocal_lattice_modes] + self.dims, 
                       dtype=complex)

        Gaussian_filter_f = self.calc_Gaussian_filter_f()

        if self.dim == 2:
                for n in range(self.number_of_primary_reciprocal_lattice_modes):
                    eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(self.psi*np.exp(
                        -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y)))

        elif self.dim == 3:
            for n in range(self.number_of_primary_reciprocal_lattice_modes):
                eta[n] = sp.fft.ifftn(Gaussian_filter_f*sp.fft.fftn(self.psi*np.exp(
                    -1j*self.q[n][0]*self.x - 1j*self.q[n][1]*self.y - 1j*self.q[n][2]*self.z  
                    )))
                
        return eta

    def calc_stress_tensor(self):
        """
        Calculates the stress of the phase-field crystal.

        Input:
            None

        Output:
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
        """ 
        Calculates the structure tensor of the phase-field crystal.

        Input:
            None

        Output:
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
        """
        Calculates the strain of the phase-field crystal.

        Input:
            None

        Output:
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

    # PLOTTING FUNCTIONS
    def plot_dislocation_nodes(self, dislocation_nodes, **kwargs):
        """
        Plots the dislocation nodes.

        Input:
            - dislocation_nodes: The dislocation nodes to plot.
            -**kwargs: Keyword arguments for the plot.
                See github.com/vidarsko/ComFiT/blob/main/docs/ClassBaseSystem.md 
                for a full list of keyword arguments.  

        Output: 
            matplotlib.axes.Axes: The axes containing the plot.
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


class PhaseFieldCrystal1DPeriodic(PhaseFieldCrystal):
    def __init__(self, nx, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal1DPeriodic'
        self.dim = 1

        # Default simulation parameters
        self.micro_resolution = kwargs.get('micro_resolution',[5])
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
        self.number_of_reciprocal_lattice_modes = 1
        self.number_of_primary_reciprocal_lattice_modes = 1

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]

        self.A = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A])

        # Set the elastic constants
        self.el_mu = 3 * self.A ** 2
        self.el_lambda = 3 * self.A ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, 
                         dx=self.dx,  dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = np.sqrt(-self.r/3 - psi0**2)
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
        self.number_of_reciprocal_lattice_modes = 3
        self.number_of_primary_reciprocal_lattice_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = np.sqrt(3) * a0 / self.micro_resolution[1]

        self.A = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A, self.A, self.A])

        # Set the elastic constants
        self.el_lambda = 3 * self.A ** 2
        self.el_mu = 3 * self.A ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)


        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        

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
    
    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((3,self.xRes,self.yRes))
        
        Lpsi = np.real(sp.fft.ifftn((1-self.calc_k2())*self.psi_f))
        stress[0] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[1]*self.dif[1]*self.psi_f))

        return stress

    
    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*self.dif[i]*field_f)*sp.fft.ifftn(self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*field_f)*sp.fft.ifftn(self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )



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
        self.number_of_reciprocal_lattice_modes = 4
        self.number_of_primary_reciprocal_lattice_modes = 2

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]

        self.A, self.B = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A, self.A, self.B, self.B])

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2
        self.el_mu = 16 * self.B ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)


        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        
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

    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((3,self.xRes,self.yRes))
        
        k2 = self.calc_k2()
        L1L2psi = np.real(sp.fft.ifftn((1-k2)*(2-k2)*self.psi_f))
        L1plusL2_f = 2 - 3*k2
        stress[0] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[1]*self.dif[1]*self.psi_f))

        return stress

    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*(2-k2)*self.dif[i]*field_f)*sp.fft.ifftn((3-2*k2)*self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*(2-k2)*field_f)*sp.fft.ifftn((3-2*k2)*self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )
            #TODO: Needs double checking (Vidar 27.02.24)


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
        self.number_of_reciprocal_lattice_modes = 6
        self.number_of_primary_reciprocal_lattice_modes = 6

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.A = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A, self.A, self.A, self.A, self.A, self.A])

        # Set the elastic constants
        self.el_lambda = 4 * self.A ** 2
        self.el_mu = 4 * self.A ** 2
        self.el_gamma = - 4*self.A**2
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True


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

    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((6,self.xRes,self.yRes,self.zRes))
        
        L1psi = np.real(sp.fft.ifftn((1-self.calc_k2())*self.psi_f))
        stress[0] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[2]*self.psi_f))
        stress[3] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[1]*self.dif[1]*self.psi_f))
        stress[4] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[1]*self.dif[2]*self.psi_f))
        stress[5] = -2*L1psi*np.real(sp.fft.ifftn(self.dif[2]*self.dif[2]*self.psi_f))

        return stress
    
    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*self.dif[i]*field_f)*sp.fft.ifftn(self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*field_f)*sp.fft.ifftn(self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )
                #TODO: Needs double checking (Vidar 27.02.24)

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
        self.number_of_reciprocal_lattice_modes = 7
        self.number_of_primary_reciprocal_lattice_modes = 4

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.A, self.B = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A, self.A, self.A, self.A, self.B, self.B, self.B])

        # Set the elastic constants
        self.el_lambda = 32/81 * self.A ** 2
        self.el_mu = 32/81 * self.A ** 2
        self.el_gamma = 32/81 * (2*self.B**2 - self.A**2)
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        
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

    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((6,self.xRes,self.yRes,self.zRes))
        
        k2 = self.calc_k2()
        L1L43psi = np.real(sp.fft.ifftn((1-k2)*(4/3-k2)*self.psi_f))
        L1plusL43_f = 7/3 - 2*k2
        stress[0] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[2]*self.psi_f))
        stress[3] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[1]*self.dif[1]*self.psi_f))
        stress[4] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[1]*self.dif[2]*self.psi_f))
        stress[5] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[2]*self.dif[2]*self.psi_f))

        return stress

    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*(4/3-k2)*self.dif[i]*field_f)*sp.fft.ifftn((7/3 - 2*k2)*self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*(4/3-k2)*field_f)*sp.fft.ifftn(7/3 - 2*k2)*self.dif[j]*(-k2)*field_f) for j in range(self.dim)])
                #TODO: Needs double checking (Vidar 27.02.24)


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
        self.number_of_reciprocal_lattice_modes = 13
        self.number_of_primary_reciprocal_lattice_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.A, self.B, self.C = self.calc_initial_amplitudes()
        self.eta0 = np.array([self.A, self.A, self.A,
                     self.B, self.B, self.B, self.B, self.B, self.B,
                     self.C, self.C, self.C, self.C])

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_mu = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_gamma = 32*self.A**2 - 16*self.B**2 - 256*self.C**2
        self.el_nu = self.el_lambda / ((self.dim - 1) * self.el_lambda + 2 * self.el_mu + self.el_gamma)

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

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

    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((6,self.xRes,self.yRes,self.zRes))
        
        k2 = self.calc_k2()
        L1L2L3psi = np.real(sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*self.psi_f))
        L2L3_p_L1L3_p_L1L2_f = 11 - 12*k2 + 3*k2**2
        stress[0] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[2]*self.psi_f))
        stress[3] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[1]*self.dif[1]*self.psi_f))
        stress[4] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[1]*self.dif[2]*self.psi_f))
        stress[5] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[2]*self.dif[2]*self.psi_f))

        return stress

    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*self.dif[i]*field_f)*sp.fft.ifftn((11-12*k2+3*k2**2)*self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*field_f)*sp.fft.ifftn((11-12*k2+3*k2**2)*self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )
        #TODO: Needs double checking (Vidar 27.02.24)