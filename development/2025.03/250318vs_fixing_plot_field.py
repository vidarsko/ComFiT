
import comfit as cf
import numpy as np
import scipy as sp

class WaveEquation(cf.BaseSystem):
    def __init__(self, dim, wave_speed, **kwargs):
        """Initialize the WaveEquation system with wave speed c."""
        self.c = wave_speed  # Wave speed
        super().__init__(dim, **kwargs)

        self.psi = np.zeros[[2] + self.dims)
        self.psi_f = np.zeros([2] + self.dims, dtype=complex)

    def conf_initial_state(self, psi):
        self.psi[0] = psi[0]
        self.psi[1] = psi[1]
        self.psi_f = self.fft(self.psi)

    def calc_omega_f(self):
        """Calculate the dispersion relation for the wave equation."""
        k2 = self.calc_k2()  # k^2 from ComFiT (sum of k[i]^2)
        return np.zeros([2] + self.dims)
    
    def calc_nonlinear_evolution_function_f(self, field, t):
        N0_f = self.fft(self.psi[1])
        N1_f = -self.c**2*self.calc_k2()*self.psi_f[0]
        return np.array([N0_f, N1_f])

    def evolve(self, number_steps):
        omega_f = self.calc_omega_f()

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, 'ETD4RK')

        for n in range(number_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, 
                                        self.calc_nonlinear_evolution_function_f, 
                                        self.psi, self.psi_f)
            self.psi = np.real(self.psi)


we = WaveEquation(1, 1)
we.conf_initial_state([we.calc_Gaussian(width=10, top=1), np.zeros(we.dims)])


fig,ax = we.plot_field(we.psi[0], ylim=[0,1])
we.show(fig)