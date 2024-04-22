import numpy as np
from comfit.core.base_system import BaseSystem
import matplotlib.pyplot as plt
from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered
import scipy as sp

class QuantumMechanics(BaseSystem):
    def __init__(self, dimension, **kwargs):
        """
        Initializes a quamtum mechanics system evolving according to the Schrödinger equation

        Input:
        - dimension : int
            The dimension of the system.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters. Same as BaseSystem

        Output:
        - QuantumMechanics object
            The system object representing the QuantumMechanics simulation.

        Example:
        qm = cf.QuantumMechanics(3,xRes=101,yRes=101,zRes=101)
        Creates a BoseEinsteinCondensate system with 3 dimensions and a spatial resolution of 101.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'QuantumMechanics1' #Quantum mechanical system with one particle.

        # Default simulation parameters
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        """
        Output a string representation of the system
        """
        return 'QuantumMechanics'

    def conf_initial_condition_Gaussian(self,position=None,width=None, initial_velocity=None):
        """
        Set the initial condition to a Gaussian wavepacket.
        """

        if position == None:
            position = self.rmid
        if width == None:
            width = self.xmax/10
    

        self.psi = np.sqrt(self.calc_Gaussian(position=position,width=width))
        
        if initial_velocity != None:
            if self.dim == 1:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * v0 * self.x)

            elif self.dim == 2:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1]*self.y))

            elif self.dim == 3:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1]*self.y + v0[2]*self.z))
        
        self.psi_f = sp.fft.fftn(self.psi)


    def conf_harmonic_potential(self,trapping_strength=None):
        """
        Set the external potential to a harmonic trap with R_tf being the thomas fermi radius
        :param R_tf:
        :Output:
        """

        if trapping_strength == None:
            trapping_strength = 1/(0.25*self.xmax)**2


        if self.dim ==1:
            return  trapping_strength*(self.x - self.xmid )**2

        if self.dim == 2:
            return trapping_strength*(((self.x-self.xmid)**2) +((self.y-self.ymid)**2) )
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2)
                                           + ((self.y - self.ymid) ** 2)
                                           +((self.z - self.zmid) ** 2) )

    def conf_hydrogen_state(self,n,l,m):
        """
        Set the initial condition to a hydrogen state with quantum numbers n,l,m
        """

        self.psi = self.calc_hydrogen_state(n,l,m)
        self.psi_f = sp.fft.fftn(self.psi)

    def calc_hydrogen_state(self,n,l,m):
        """
        Calculate the hydrogen state with quantum numbers n,l,m

        Input:
            n: principal quantum number
            l: azimuthal quantum number
            m: magnetic quantum number

        Output:
            (np.ndarray): wavefunction
        """

        if not self.dim == 3:
            raise(Exception('The hydrogen state is only defined in 3D'))

        r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        theta = np.arccos(self.z/r)
        theta[np.isnan(theta)] = 0
        phi = np.arctan2(self.y,self.x)

        R = (2*r/n)**l * np.exp(-r/n) *  sp.special.genlaguerre(n-l-1,2*l+1)(2*r/n)
        Y = sp.special.sph_harm(m,l,phi,theta)
        
        return R*Y

    

    ## Calculation functions
    def calc_nonlinear_evolution_function_f(self,psi,t):

        return sp.fft.fftn((1j) * (-self.V_ext) * psi)


    ## Time evolution functions
    def evolve_schrodinger(self,number_of_steps,method = 'ETD2RK'):
        omega_f = -1j / 2 * self.calc_k2()
        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop
        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise(Exception('This method is not implemented'))
        for n in range(number_of_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, self.calc_nonlinear_evolution_function_f,
                                                           self.psi, self.psi_f)


    # Hamilton minimization functions

    def calc_Hamiltonian(self):
        integrand = -1/2*np.conj(self.psi)\
            * sp.fft.ifftn(-self.calc_k2()*self.psi_f)\
            + self.V_ext*abs(self.psi)**2
        return self.calc_integrate_field(integrand)
    
    ## Plotting functions 
    def plot(self, ylim=None):

        if self.dim == 1:
            plt.clf()
            ax = plt.gcf().add_subplot(111)

            cmap = tool_colormap_angle()


            y = np.abs(self.psi) ** 2

            c = np.angle(self.psi)

            ax.plot(self.x,y,color='black')

            if ylim is not None:
                ax.set_ylim(0,ylim)


            # Color in the graph based on the argument of the wavefunction
            alpha_level=0.7

            ax.fill_between([self.xmin,self.xmin+self.dx/2], [y[0],(y[0]+y[1])/2],
                            color=cmap((c[0] + np.pi) / (2 * np.pi)), alpha=alpha_level,edgecolor='none')

            for i in range(1,self.xRes-1):
                ax.fill_between([self.x[i]-self.dx/2,self.x[i]+self.dx/2], [(y[i]+y[i-1])/2,(y[i]+y[i+1])/2],
                                color = cmap((c[i] + np.pi) / (2 * np.pi)), alpha=alpha_level,edgecolor='none')

            ax.fill_between([self.xmax-self.dx/2,self.xmax], [y[-1],(y[-1]+y[-2])/2],
                            color=cmap((c[-1] + np.pi) / (2 * np.pi)), alpha=alpha_level,edgecolor='none')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$|\psi|^2$')

            sm = plt.cm.ScalarMappable(cmap=cmap)
            cbar = plt.gcf().colorbar(sm, ax=ax)

            #cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
            cbar.set_ticks(np.array([0, 1/6, 2/6, 3/6, 4/6, 5/6, 1]))
            cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])
            cbar.set_label('arg($\psi$)')


            return ax
    