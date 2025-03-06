import comfit as cf
import numpy as np
import scipy as sp

class TwoParticleQM(cf.BaseSystem):
    def __init__(self, **kwargs):
        """Initialize a 2-particle QM system.
        
        Args:
            **kwargs: Additional arguments for BaseSystem (e.g., xRes, dx, dt)
        """
        super().__init__(2, **kwargs)
        
        # Initialize wavefunction (psi) and its Fourier transform
        self.psi = np.zeros(2, dtype=complex)  # Shape: (xRes, yRes) for 2D
        self.psi_f = np.zeros(2, dtype=complex)
        
        # Coordinate arrays (x1 = x, x2 = y in ComFiT convention)
        self.x1 = self.x  # Shape: (xRes, 1)
        self.x2 = self.y  # Shape: (1, yRes)

    def conf_initial_condition(self, psi_initial):
        """Set initial wavefunction."""
        self.psi = psi_initial
        self.psi_f = self.fft(self.psi)


    def calc_omega_f(self):
        """Calculate linear evolution operator in Fourier space.
        H = -ħ²/(2m1)∇₁² - ħ²/(2m2)∇₂² + V(x1, x2)
        Fourier transform of kinetic terms: -ħ²/(2m)k²
        """
        # k[0] is wave numbers for x1, k[1] for x2
        k1_sq = self.k[0]**2  # Shape: (xRes, 1)
        k2_sq = self.k[1]**2  # Shape: (1, yRes)
        
        # Kinetic terms in Fourier space
        return  1j/2 * (k1_sq +  k2_sq)
    
    def calc_nonlinear_evolution_function_f(self, field, t):
        return self.fft(self.V_ext*self.psi)

    def evolve(self, number_steps):
        omega_f = self.calc_omega_f()

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, 'ETD2RK')

        for n in range(number_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, 
                                        self.calc_nonlinear_evolution_function_f, 
                                        self.psi, self.psi_f)

    def conf_initial_condition_Gaussian(self,
            particle_type,
            position1,
            width1,
            initial_velocity1,
            position2,
            width2,
            initial_velocity2):

        psi1 = np.sqrt(self.calc_Gaussian(self.x1, position1, width1))*np.exp(1j * initial_velocity1 * self.x1)  \
              * np.sqrt(self.calc_Gaussian(self.x2, position2, width2))*np.exp(1j * initial_velocity2 * self.x2)

        psi2 = np.sqrt(self.calc_Gaussian(self.x2, position1, width1))*np.exp(1j * initial_velocity1 * self.x2)  \
              * np.sqrt(self.calc_Gaussian(self.x1, position2, width2))*np.exp(1j * initial_velocity2 * self.x1)


        if particle_type == 'boson':
            self.psi = 1/np.sqrt(2)*(psi1 + psi2)
        elif particle_type == 'fermion':
            self.psi = 1/np.sqrt(2)*(psi1 - psi2)
        
        self.psi_f = sp.fft.fftn(self.psi)

    def calc_Gaussian(self, x, position, width):

        rx2m = (x - position - self.size_x) ** 2
        rx2 = (x - position) ** 2
        rx2p = (x - position + self.size_x) ** 2

        r2 = np.min(np.stack((rx2m, rx2, rx2p)), axis=0)

        return (2*np.pi*width**2)**(-1/2)*np.exp(-r2/(2*width**2))

if __name__ == '__main__':
    # Initialize system
    qm = TwoParticleQM(xlim=[-20,20], ylim=[-20,20])
    qm.conf_initial_condition_Gaussian(particle_type='fermion',
                                           position1=-10, width1=3, initial_velocity1=-20,
                                           position2=10, width2=3, initial_velocity2=0)
    

    qm.V_ext = -qm.calc_Gaussian(qm.x1-qm.x2, 0, 2)
    # fig, ax = qm.plot_complex_field(qm.psi)
    # qm.show(fig)


    qm1 = cf.BaseSystem(1, xlim=qm.xlim)

    psi_proj = np.sum(abs(qm.psi)**2, axis=0)
    ymax = np.max(psi_proj)

    for n in range(100):
        qm.evolve(1)
        psi_proj = np.sum(abs(qm.psi)**2, axis=0)

        # fig, ax = qm.plot_complex_field(qm.psi)
        fig, ax = qm1.plot_field(psi_proj, ylim=[0, ymax])
        qm.plot_save(fig,n)
    cf.tool_make_animation_gif(n)

    # print(qm.calc_integrate_field(abs(qm.psi)**2))
