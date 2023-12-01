import numpy as np
from comfit.core.base_system import BaseSystem
import matplotlib.pyplot as plt
from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered

class QM(BaseSystem):
    def __init__(self, dimension, **kwargs):
        """
        Initializes a quamtum mechanics system evolving according to the Schrödinger equation

        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'QM1' #Quantum mechanical system with one particle.

        # Default simulation parameters
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    def set_initial_condition_gaussian(self,position=None,width=None):

        if self.dim == 1:
            if position == None:
                position = self.xmid
            if width == None:
                width = 5

            self.psi = np.sqrt(1/np.sqrt(2*np.pi*width**2)*np.exp(-(self.x-position)**2/(2*width**2)))
            self.psi_f = np.fft.fftn(self.psi)



    def set_harmonic_potential(self,trapping_strength=None):
        """
        Set the external potential to a harmonic trap with R_tf being the thomas fermi radius
        :param R_tf:
        :return:
        """

        if trapping_strength == None:
            trapping_strength = 1/(0.25*self.xmax)**2


        if self.dim ==1:
            return  trapping_strength*(self.x - self.xmid )**2

        if self.dim == 2:
            return trapping_strength*(((self.x-self.xmid)**2).reshape(self.xRes, 1)
                                         +((self.y-self.ymid)**2).reshape(1, self.yRes) )
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1,1)
                                           + ((self.y - self.ymid) ** 2).reshape(1, self.yRes,1)
                                           +((self.z - self.zmid) ** 2).reshape(1, 1,self.zRes) )

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



    def calc_nonlinear_evolution_function_f(self,psi):

        return np.fft.fftn((1j) * (-self.V_ext) * psi)


    def evolve_schrodinger(self,number_of_steps):

        omega_f = -1j/2*self.calc_k2()

        integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)

        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETD2RK_loop(integrating_factors_f, self.calc_nonlinear_evolution_function_f,
                                                           self.psi, self.psi_f)


    # Hamilton minimization functions

    def calc_Hamiltonian(self):
        integrand = -1/2*np.conj(self.psi)\
            * np.fft.ifftn(-self.calc_k2()*self.psi_f)\
            + self.V_ext*abs(self.psi)**2
        return self.calc_integrate_field(integrand)
    