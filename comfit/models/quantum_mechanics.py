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

            self.psi = np.sqrt(1/np.sqrt(2*np.pi*width**2)*np.exp(-(self.x-self.xmid)**2/(2*width**2)))
            self.psi_f = np.fft.fftn(self.psi)

    def plot(self):

        self.plot_field(np.abs(self.psi)**2,1) #This ax argument is just a dummy atm. Needs to be fixed (Vidar 14.09.23)



    def calc_nonlinear_evolution_term_f(self):


    def evolve_schrodinger(self):
