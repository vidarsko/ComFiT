import numpy as np
from comfit.core.base_system import BaseSystem

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
            return  trapping_strength*(self.x -self.xmid )**2
        if self.dim == 2:
            return trapping_strength*(((self.x-self.xmid)**2).reshape(self.xRes, 1)
                                         +((self.y-self.ymid)**2).reshape(1, self.yRes) )
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1,1)
                                           + ((self.y - self.ymid) ** 2).reshape(1, self.yRes,1)
                                           +((self.z - self.zmid) ** 2).reshape(1, 1,self.zRes) )

    def plot(self):

        self.plot_field(np.abs(self.psi)**2,1) #This ax argument is just a dummy atm. Needs to be fixed (Vidar 14.09.23)



    def calc_nonlinear_evolution_term_f(self,psi):

        return np.fft.fftn((1j) * (-self.V_ext) * psi)


    def calc_evolution_integrating_factors_schrodinger_f(self):

        k2 = self.calc_k2()

        omega_f = -1j *  1 / 2 * k2

        integration_factors_f = [0,0,0]

        integration_factors_f[0] = np.exp(omega_f * self.dt)
        If1 = integration_factors_f[0]


        integration_factors_f[1] = (If1 - 1) / omega_f

        integration_factors_f[2] = 1 / (self.dt * omega_f**2) * (If1 - 1 - omega_f * self.dt)

        for i in range(3):
            if self.dim == 1:
                integration_factors_f[i][0]=0
            elif self.dim == 2:
                integration_factors_f[i][0,0]=0
            elif self.dim == 3:
                integration_factors_f[i][0,0,0]=0

        return integration_factors_f

    def evolve_schrodinger(self,number_of_steps):

        integrating_factors_f = self.calc_evolution_integrating_factors_schrodinger_f()

        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK2_loop(integrating_factors_f, self.psi, self.psi_f)
