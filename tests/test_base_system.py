import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

class TestBaseSystem(unittest.TestCase):

    def test_init_empty(self):
        """ Test BaseSystem initialization with  parameters """
        for dim in [1,2,3]:
            try:
                bs = cf.BaseSystem(dim)
                self.assertIsInstance(bs, cf.BaseSystem)
            except Exception as e:
                self.fail(f"Initialization failed with no parameters: {e}")

    def test_init_with_dim_1_and_params(self):
        """ Test 1D BaseSystem initialization with different parameters """
        params = [
            {'xRes': 11, 'dx': 0.1},
            {'xRes': 32, 'dx': 0.3},
            # Add more parameter combinations as needed
        ]
        for p in params:
            with self.subTest(p=p):
                try:
                    bs = cf.BaseSystem(1,**p)
                    self.assertIsInstance(bs, cf.BaseSystem)
                except Exception as e:
                    self.fail(f"Initialization failed with parameters {p}: {e}")

    def test_init_with_dim_2_and_params(self):
        """ Test 2D BaseSystem initialization with different parameters """
        params = [
            {'xRes': 11, 'dx': 0.1, 'yRes': 16, 'dy': 0.1},
            {'xRes': 28, 'dx': 0.2, 'yRes': 32, 'dy': 1.3},
            # Add more parameter combinations as needed
        ]
        for p in params:
            with self.subTest(p=p):
                try:
                    bs = cf.BaseSystem(2,**p)
                    self.assertIsInstance(bs, cf.BaseSystem)
                except Exception as e:
                    self.fail(f"Initialization failed with parameters {p}: {e}")

    def test_init_with_dim_3_and_params(self):
        """ Test 3D BaseSystem initialization with different parameters """
        params = [
            {'xRes': 11, 'dx': 0.1, 'yRes': 16, 'dy': 0.1, 'zRes': 16, 'dz': 0.1},
            {'xRes': 28, 'dx': 0.2, 'yRes': 32, 'dy': 1.3, 'zRes': 32, 'dz': 0.2},
            # Add more parameter combinations as needed
        ]
        for p in params:
            with self.subTest(p=p):
                try:
                    bs = cf.BaseSystem(3,**p)
                    self.assertIsInstance(bs, cf.BaseSystem)
                except Exception as e:
                    self.fail(f"Initialization failed with parameters {p}: {e}")

    def test_str(self):
        """ Test BaseSystem string representation """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1)

        # Get string representation
        result = bs.__str__()

        # Check if result is a string
        assert isinstance(result, str), "Expected string representation of BaseSystem object"

    def test_calc_angle_field_single_vortex(self):
        """ Test calculation of angle field for single vortex """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Calculate angle field
        angle_field = bs.calc_angle_field_single_vortex(position=bs.rmid,charge=-1)

        # Check if angle field is a numpy array
        assert isinstance(angle_field, np.ndarray), "Expected numpy array for angle field"

        # Check if angle field has the correct shape
        assert angle_field.shape == (bs.yRes, bs.xRes), "Expected angle field with shape (yRes, xRes)"

        # Check if angle is between -pi and pi
        assert np.all(angle_field >= -np.pi) and np.all(angle_field <= np.pi), "Expected angle field between -pi and pi"
        

    def test_calc_angle_field_vortex_dipole(self):
        """ Test calculation of angle field for vortex dipole """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Calculate angle field
        angle_field = bs.calc_angle_field_vortex_dipole(dipole_vector=[0.24*bs.xmax,0.3*bs.ymax], dipole_position=bs.rmid)

        # Check if angle field is a numpy array
        assert isinstance(angle_field, np.ndarray), "Expected numpy array for angle field"

        # Check if angle field has the correct shape
        assert angle_field.shape == (bs.yRes, bs.xRes), "Expected angle field with shape (yRes, xRes)"

        # Check if angle is between -pi and pi
        assert np.all(angle_field >= -np.pi) and np.all(angle_field <= np.pi), "Expected angle field between -pi and pi"

    def test_calc_angle_field_vortex_ring(self):
        """ Test calculation of angle field for vortex ring """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Calculate angle field
        angle_field = bs.calc_angle_field_vortex_ring(normal_vector=[0.1,0.4,0.2], position=bs.rmid)

        # Check if angle field is a numpy array
        assert isinstance(angle_field, np.ndarray), "Expected numpy array for angle field"

        # Check if angle field has the correct shape
        assert angle_field.shape == (bs.zRes, bs.yRes, bs.xRes), "Expected angle field with shape (zRes, yRes, xRes)"

        # Check if angle is between -pi and pi
        assert np.all(angle_field >= -np.pi) and np.all(angle_field <= np.pi), "Expected angle field between -pi and pi"


    def test_calc_wavenums(self):
        """ Test calculation of wavenumbers """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1)

        x = np.array([-10, -5, 0, 5, 10])
        k = bs.calc_wavenums(x)

        #Check values of k
        np.testing.assert_allclose(k, np.array([ 0., 0.25132741,0.50265482,-0.50265482,-0.25132741]))


    def test_calc_determinant_field(self):
        """ Test calculation of determinant field """
        
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=101, dx=1, yRes=101, dy=1)

        # 2 dimensions
        # Testing a zero-field
        psi = np.zeros((2,bs.xRes,bs.yRes))
        D = bs.calc_determinant_field(psi)
        np.testing.assert_allclose(D, np.zeros((bs.xRes,bs.yRes)))

        # Testing a +1 charge field
        # Making the field respecting the periodicity
        psi[0] = bs.xmax/(2*np.pi)*np.sin((bs.x - bs.xmid)/bs.xmax*2*np.pi)
        psi[1] = bs.ymax/(2*np.pi)*np.sin((bs.y - bs.ymid)/bs.ymax*2*np.pi)
        D = bs.calc_determinant_field(psi)
        self.assertAlmostEquals(D[bs.xmidi,bs.ymidi], 1.0)

        # Testing a -1 charge field
        # Making the field respecting the periodicity
        psi[0] = bs.xmax/(2*np.pi)*np.sin((bs.x - bs.xmid)/bs.xmax*2*np.pi)
        psi[1] = -bs.ymax/(2*np.pi)*np.sin((bs.y - bs.ymid)/bs.ymax*2*np.pi)
        D = bs.calc_determinant_field(psi)
        self.assertAlmostEquals(D[bs.xmidi,bs.ymidi], -1.0)

        # 3 dimensions
        # Testing a 0-field
        bs = cf.BaseSystem(3, xRes=31, dx=1, yRes=42, dy=1, zRes=13, dz=1)
        psi = np.zeros((2,bs.xRes,bs.yRes,bs.zRes))
        D = bs.calc_determinant_field(psi)
        np.testing.assert_allclose(D, np.zeros((3,bs.xRes,bs.yRes,bs.zRes)))
    
        # Testing a +1 charge field
        # Making the field respecting the periodicity
        psi[0] = bs.xmax/(2*np.pi)*np.sin((bs.x - bs.xmid)/bs.xmax*2*np.pi)
        psi[1] = bs.zmax/(2*np.pi)*np.sin((bs.z - bs.zmid)/bs.zmax*2*np.pi)
        D = bs.calc_determinant_field(psi)
        # At the center, the defect density should point in the negative y-direction
        np.testing.assert_allclose(np.array([
                                 D[0,bs.xmidi,bs.ymidi,bs.zmidi],
                                 D[1,bs.xmidi,bs.ymidi,bs.zmidi],
                                 D[2,bs.xmidi,bs.ymidi,bs.zmidi]]), 
                                 np.array([0.0, -1.0, 0.0]),atol=1e-15)



    def test_evolution_scheme_1d(self):
        """ Test 1D evolution scheme """
        
        # Testing the wave equation in one dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=101, dx=1)

        # Set test parameters
        bs.A = 0.1
        bs.sigma = 5
        bs.T0 = 20

        # Set end time
        end_time = 200

        # Scipy benchmark

        # Initial condition
        T_initial = np.zeros((bs.xRes))

        # Forcing term function
        def forcing_term(T):
            return bs.A*(bs.T0-T)*np.exp(-(bs.x-bs.xmid)**2/(2*bs.sigma**2))
        
        # Function representing the discretized PDE
        def heat_equation(t, T):
            dTdx2 = np.roll(T, -1) - 2 * T + np.roll(T, 1)
            dTdx2 /= bs.dx**2
            return dTdx2 + forcing_term(T)

        # Time span for the integration
        t_span = (0, end_time)

        # Solve the equation
        T_benchmark = sp.integrate.solve_ivp(heat_equation, t_span, T_initial, method='RK45')

        # ComFiT simulation

        for method in ['ETD2RK','ETD4RK']:
            # ComFiT simulation

            # Add linear and non-linear functions
            def calc_omega_f(bs):
                return -bs.calc_k2()

            def calc_nonlinear_evolution_function_f(bs,T,t):
                f = bs.A*(bs.T0-T)*np.exp(-(bs.x-bs.xmid)**2/(2*bs.sigma**2))
                return sp.fft.fft(f)

            # Add functions to BaseSystem object
            bs.calc_omega_f = calc_omega_f.__get__(bs)
            bs.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f.__get__(bs)

            # Define evolve function
            def evolve(bs,number_of_steps, method = method):
                omega_f = bs.calc_omega_f()

                integrating_factors_f, solver = bs.calc_integrating_factors_f_and_solver(omega_f,method)

                for n in range(number_of_steps):
                    bs.T, bs.T_f = solver(integrating_factors_f,
                                            bs.calc_nonlinear_evolution_function_f,
                                            bs.T, bs.T_f)
                    bs.T = np.real(bs.T)
            
            # Add evolve function to BaseSystem object
            bs.evolve = evolve.__get__(bs)

            # Initial condition
            bs.T = np.zeros((bs.xRes))
            bs.T_f = sp.fft.fft(bs.T)   

            # Evolve the system
            number_of_steps = int(end_time/bs.dt)
            bs.evolve(number_of_steps)

            # Check if the solution is close to the benchmark
            np.testing.assert_allclose(bs.T, T_benchmark.y[:,-1], atol=1e-2, rtol=1e-2)
            
        
    def test_evolution_scheme_2d(self):
        """ Test 2D evolution scheme"""

        # Testing the wave equation in two dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=51, dx=1, yRes=51, dy=1)

        # Set test parameters
        bs.A = 0.1
        bs.sigma = 2
        bs.T0 = 20

        # Set end time
        end_time = 100

        # Scipy benchmark

        # Initial condition
        T_initial = np.zeros((bs.xRes,bs.yRes)).flatten()

        # Forcing term function
        def forcing_term(T_flat):
            T = T_flat.reshape((bs.xRes, bs.yRes))
            forcing = bs.A * (bs.T0 - T) * np.exp(-((bs.x - bs.xmid) ** 2 + (bs.y - bs.ymid) ** 2) / (2 * bs.sigma ** 2))
            return forcing.flatten()

        # Function representing the discretized PDE
        def heat_equation(t,T_flat):
            T = T_flat.reshape((bs.xRes, bs.yRes))
            dTdx2 = np.roll(T, -1, axis=0) - 2 * T + np.roll(T, 1, axis=0)
            dTdx2 /= bs.dx ** 2
            dTdy2 = np.roll(T, -1, axis=1) - 2 * T + np.roll(T, 1, axis=1)
            dTdy2 /= bs.dy ** 2
            return dTdx2.flatten() + dTdy2.flatten() + forcing_term(T_flat)
        
        # Time span for the integration
        t_span = (0, end_time)

        # Solve the equation
        T_benchmark = sp.integrate.solve_ivp(heat_equation, t_span, T_initial, method='RK45')
        T_final_2D = T_benchmark.y[:, -1].reshape((bs.xRes, bs.yRes))

        # ComFiT simulation

        for method in ['ETD2RK','ETD4RK']:

            # Add linear and non-linear functions
            def calc_omega_f(bs):
                return -bs.calc_k2()

            def calc_nonlinear_evolution_function_f(bs,T,t):
                f = bs.A*(bs.T0-T)*np.exp(-((bs.x-bs.xmid)**2+(bs.y-bs.ymid)**2)/(2*bs.sigma**2))
                return sp.fft.fft2(f)
            
            # Add functions to BaseSystem object
            bs.calc_omega_f = calc_omega_f.__get__(bs)
            bs.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f.__get__(bs)

            # Define evolve function
            def evolve(bs,number_of_steps, method = method):
                omega_f = bs.calc_omega_f()

                integrating_factors_f, solver = bs.calc_integrating_factors_f_and_solver(omega_f,method)

                for n in range(number_of_steps):
                    bs.T, bs.T_f = solver(integrating_factors_f,
                                            bs.calc_nonlinear_evolution_function_f,
                                            bs.T, bs.T_f)
                    bs.T = np.real(bs.T)
            
            # Add evolve function to BaseSystem object
            bs.evolve = evolve.__get__(bs)

            # Initial condition
            bs.T = np.zeros((bs.xRes,bs.yRes))
            bs.T_f = sp.fft.fft2(bs.T)

            # Evolve the system
            number_of_steps = int(end_time/bs.dt)
            bs.evolve(number_of_steps)

            # Check if the solution is close to the benchmark
            np.testing.assert_allclose(bs.T, T_final_2D, atol=1e-2, rtol=1e-2)

    def test_evolution_scheme_3d(self):
        """ Test 3D evolution scheme """

        # Testing the wave equation in three dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=21, dx=1, yRes=21, dy=1, zRes=21, dz=1)

        # Set test parameters
        bs.A = 0.1
        bs.sigma = 2
        bs.T0 = 20

        # Set end time
        end_time = 100

        # Scipy benchmark

        # Initial condition
        T_initial = np.zeros((bs.xRes,bs.yRes,bs.zRes)).flatten()

        # Forcing term function
        def forcing_term(T_flat):
            T = T_flat.reshape((bs.xRes, bs.yRes, bs.zRes))
            forcing = bs.A * (bs.T0 - T) * np.exp(-((bs.x - bs.xmid) ** 2 + (bs.y - bs.ymid) ** 2 + (bs.z - bs.zmid) ** 2) / (2 * bs.sigma ** 2))
            return forcing.flatten()
        
        # Function representing the discretized PDE
        def heat_equation(t,T_flat):
            T = T_flat.reshape((bs.xRes, bs.yRes, bs.zRes))
            dTdx2 = np.roll(T, -1, axis=0) - 2 * T + np.roll(T, 1, axis=0)
            dTdx2 /= bs.dx ** 2
            dTdy2 = np.roll(T, -1, axis=1) - 2 * T + np.roll(T, 1, axis=1)
            dTdy2 /= bs.dy ** 2
            dTdz2 = np.roll(T, -1, axis=2) - 2 * T + np.roll(T, 1, axis=2)
            dTdz2 /= bs.dz ** 2
            return dTdx2.flatten() + dTdy2.flatten() + dTdz2.flatten() + forcing_term(T_flat)
        
        # Time span for the integration
        t_span = (0, end_time)

        # Solve the equation
        T_benchmark = sp.integrate.solve_ivp(heat_equation, t_span, T_initial, method='RK45')
        T_final_3D = T_benchmark.y[:, -1].reshape((bs.xRes, bs.yRes, bs.zRes))

        # ComFiT simulation

        for method in ['ETD2RK','ETD4RK']:

            # Add linear and non-linear functions
            def calc_omega_f(bs):
                return -bs.calc_k2()

            def calc_nonlinear_evolution_function_f(bs,T,t):
                f = bs.A*(bs.T0-T)*np.exp(-((bs.x-bs.xmid)**2+(bs.y-bs.ymid)**2+(bs.z-bs.zmid)**2)/(2*bs.sigma**2))
                return sp.fft.fftn(f)
            
            # Add functions to BaseSystem object
            bs.calc_omega_f = calc_omega_f.__get__(bs)
            bs.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f.__get__(bs)

            # Define evolve function
            def evolve(bs,number_of_steps, method = 'ETD2RK'):
                omega_f = bs.calc_omega_f()

                integrating_factors_f, solver = bs.calc_integrating_factors_f_and_solver(omega_f,method)

                for n in range(number_of_steps):
                    bs.T, bs.T_f = solver(integrating_factors_f,
                                            bs.calc_nonlinear_evolution_function_f,
                                            bs.T, bs.T_f)
                    bs.T = np.real(bs.T)

            # Add evolve function to BaseSystem object
            bs.evolve = evolve.__get__(bs)

            # Initial condition
            bs.T = np.zeros((bs.xRes,bs.yRes,bs.zRes))
            bs.T_f = sp.fft.fftn(bs.T)

            # Evolve the system
            number_of_steps = int(end_time/bs.dt)
            bs.evolve(number_of_steps)

            # Check if the solution is close to the benchmark
            np.testing.assert_allclose(bs.T, T_final_3D, atol=1e-2, rtol=1e-2)

        
    def test_plot_field(self):
        """ Test plotting of field """

        # 1 dimension
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=32, dx=0.1)

        # Create field
        field = np.random.rand(bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            plt.close()
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 2 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Create field
        field = np.random.rand(bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            plt.close()
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 3 dimensions
        # Initialize BaseSystem object

        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Create field
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            plt.close()
        except Exception as e:
            self.fail(f"Plotting failed: {e}")
    


        


if __name__ == '__main__':
    # unittest.main()

    unittest.main(defaultTest='TestBaseSystem.test_plot_field', verbosity=2)
