import unittest
import sys
import os

# Adjust the path to import your comfit package
sys.path.append(os.path.abspath('../'))

import comfit as cf

class TestBaseSystem(unittest.TestCase):

    def test_init_empty(self):
        """ Test BaseSystem initialization with  parameters """
        for dim in range(3):
            try:
                bs = cf.BaseSystem(dim)
                self.assertIsInstance(bs, cf.BaseSystem)
            except Exception as e:
                self.fail(f"Initialization failed with no parameters: {e}")

    # def test_init_with_params(self):
    #     """ Test BaseSystem initialization with different parameters """
    #     params = [
    #         {'param1': 'value1', 'param2': 10},
    #         {'param1': 'value2', 'param2': 20},
    #         # Add more parameter combinations as needed
    #     ]
    #     for p in params:
    #         with self.subTest(p=p):
    #             try:
    #                 bs = cf.BaseSystem(**p)
    #                 self.assertIsInstance(bs, cf.BaseSystem)
    #             except Exception as e:
    #                 self.fail(f"Initialization failed with parameters {p}: {e}")

if __name__ == '__main__':
    unittest.main()
