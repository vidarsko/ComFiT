import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np

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


    # def test_calc_wavenums(self):
    #     """ Test calculation of wavenumbers """
        
    #     # Initialize BaseSystem object
    #     bs = cf.BaseSystem(2)

    #     x = np.array([-10, -5, 0, 5, 10])
    #     k = bs.calc_wavenums(x)

    #     #Check values of k

        
if __name__ == '__main__':
    unittest.main()
