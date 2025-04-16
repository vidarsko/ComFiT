import unittest
import numpy as np
import comfit as cf  # Assuming comfit and the tool_configure_axis function are available

from comfit.tool import tool_configure_axis as cf_tool_configure_axis

# Default values used within the tool_configure_axis function for testing comparison
DEFAULT_XRES = 101
DEFAULT_DX = 1.0
DEFAULT_XMIN = 0
DEFAULT_XMAX = 101

class TestConfigureAxis(unittest.TestCase):

    def assertAxisConfig(self, dim, name, args, expected_xlim, expected_xmin, expected_xmax, expected_xRes, expected_dx):
        """Helper method to run a test case and assert results."""
        xlim, xmin, xmax, xRes, dx = cf.tool.tool_configure_axis(dim, name, **args)
        np.testing.assert_allclose(xlim, expected_xlim, rtol=1e-7, atol=1e-7, err_msg=f"xlim mismatch for args={args}")
        self.assertAlmostEqual(xmin, expected_xmin, places=7, msg=f"xmin mismatch for args={args}")
        self.assertAlmostEqual(xmax, expected_xmax, places=7, msg=f"xmax mismatch for args={args}")
        self.assertEqual(xRes, expected_xRes, msg=f"xRes mismatch for args={args}")
        self.assertAlmostEqual(dx, expected_dx, places=7, msg=f"dx mismatch for args={args}")

    # Test cases primarily for 'x' axis in dim=2 (y/z defaults don't interfere)
    def test_no_values_provided(self):
        self.assertAxisConfig(dim=2, name='x', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)

    def test_only_dx(self):
        dx_val = 0.1
        expected_xRes = round((DEFAULT_XMAX - DEFAULT_XMIN) / dx_val)
        self.assertAxisConfig(dim=2, name='x', args={'dx': dx_val},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - dx_val],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=expected_xRes, expected_dx=dx_val)

    def test_only_xRes(self):
        xRes_val = 56
        expected_xmax = xRes_val
        expected_dx = (expected_xmax - DEFAULT_XMIN) / xRes_val
        self.assertAxisConfig(dim=2, name='x', args={'xRes': xRes_val},
                              expected_xlim=[DEFAULT_XMIN, expected_xmax - expected_dx],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_xRes_and_dx(self):
        xRes_val = 56
        dx_val = 0.1
        expected_xmax = DEFAULT_XMIN + xRes_val * dx_val
        self.assertAxisConfig(dim=2, name='x', args={'xRes': xRes_val, 'dx': dx_val},
                              expected_xlim=[DEFAULT_XMIN, expected_xmax - dx_val],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=dx_val)

    def test_only_xmax(self):
        xmax_val = 56.0
        expected_xmin = DEFAULT_XMIN # xmax > 0
        expected_dx = (xmax_val - expected_xmin) / DEFAULT_XRES
        self.assertAxisConfig(dim=2, name='x', args={'xmax': xmax_val},
                              expected_xlim=[expected_xmin, xmax_val - expected_dx],
                              expected_xmin=expected_xmin, expected_xmax=xmax_val,
                              expected_xRes=DEFAULT_XRES, expected_dx=expected_dx)

    def test_xmax_and_dx(self):
        xmax_val = 21.0
        dx_val = 0.1
        expected_xmin = DEFAULT_XMIN # xmax > 0
        expected_xRes = round((xmax_val - expected_xmin) / dx_val)
        self.assertAxisConfig(dim=2, name='x', args={'xmax': xmax_val, 'dx': dx_val},
                              expected_xlim=[expected_xmin, xmax_val - dx_val],
                              expected_xmin=expected_xmin, expected_xmax=xmax_val,
                              expected_xRes=expected_xRes, expected_dx=dx_val)

    def test_xmax_and_xRes(self):
        xmax_val = 21.0
        xRes_val = 56
        expected_xmin = DEFAULT_XMIN # xmax > 0
        expected_dx = (xmax_val - expected_xmin) / xRes_val
        self.assertAxisConfig(dim=2, name='x', args={'xmax': xmax_val, 'xRes': xRes_val},
                              expected_xlim=[expected_xmin, xmax_val - expected_dx],
                              expected_xmin=expected_xmin, expected_xmax=xmax_val,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_xmax_xRes_dx(self):
        xmax_val = 21.0
        xRes_val = 56
        dx_val = 0.1
        expected_xmin = xmax_val - xRes_val * dx_val
        self.assertAxisConfig(dim=2, name='x', args={'xmax': xmax_val, 'xRes': xRes_val, 'dx': dx_val},
                              expected_xlim=[expected_xmin, xmax_val - dx_val],
                              expected_xmin=expected_xmin, expected_xmax=xmax_val,
                              expected_xRes=xRes_val, expected_dx=dx_val)

    def test_only_xmin(self):
        xmin_val = -5.0
        expected_xmax = xmin_val + DEFAULT_XMAX
        expected_dx = DEFAULT_DX
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val},
                              expected_xlim=[xmin_val, expected_xmax - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=expected_xmax,
                              expected_xRes=DEFAULT_XRES, expected_dx=expected_dx)

    def test_xmin_and_dx(self):
        xmin_val = -5.0
        dx_val = 0.1
        expected_xRes = DEFAULT_XRES
        expected_xmax = xmin_val + expected_xRes * dx_val
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'dx': dx_val},
                              expected_xlim=[xmin_val, expected_xmax - dx_val],
                              expected_xmin=xmin_val, expected_xmax=expected_xmax,
                              expected_xRes=expected_xRes, expected_dx=dx_val)

    def test_xmin_and_xRes(self):
        xmin_val = -5.0
        xRes_val = 56
        expected_dx = DEFAULT_DX
        expected_xmax = xmin_val + xRes_val * expected_dx
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xRes': xRes_val},
                              expected_xlim=[xmin_val, expected_xmax - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_xmin_xRes_dx(self):
        xmin_val = -5.0
        xRes_val = 56
        dx_val = 0.1
        expected_xmax = xmin_val + xRes_val * dx_val
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xRes': xRes_val, 'dx': dx_val},
                              expected_xlim=[xmin_val, expected_xmax - dx_val],
                              expected_xmin=xmin_val, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=dx_val)

    def test_xmin_and_xmax(self):
        xmin_val = -5.0
        xmax_val = 21.0
        expected_xRes = round(xmax_val - xmin_val) # Uses default dx=1 for calculation
        expected_dx = DEFAULT_DX
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xmax': xmax_val},
                              expected_xlim=[xmin_val, xmax_val - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=xmax_val,
                              expected_xRes=expected_xRes, expected_dx=expected_dx)

    def test_xmin_xmax_dx(self):
        xmin_val = -5.0
        xmax_val = 21.0
        dx_val = 0.1
        expected_xRes = round((xmax_val - xmin_val) / dx_val)
        # dx gets recalculated to fit exactly
        expected_dx = (xmax_val - xmin_val) / expected_xRes
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xmax': xmax_val, 'dx': dx_val},
                              expected_xlim=[xmin_val, xmax_val - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=xmax_val,
                              expected_xRes=expected_xRes, expected_dx=expected_dx)

    def test_xmin_xmax_xRes(self):
        xmin_val = -5.0
        xmax_val = 21.0
        xRes_val = 56
        expected_dx = (xmax_val - xmin_val) / xRes_val
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xmax': xmax_val, 'xRes': xRes_val},
                              expected_xlim=[xmin_val, xmax_val - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=xmax_val,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_xmin_xmax_xRes_dx_hierarchy(self):
        # dx should be ignored and recalculated based on xmin, xmax, xRes
        xmin_val = -5.0
        xmax_val = 21.0
        xRes_val = 56
        dx_ignored = 0.1
        expected_dx = (xmax_val - xmin_val) / xRes_val # dx recalculated
        self.assertAxisConfig(dim=2, name='x', args={'xmin': xmin_val, 'xmax': xmax_val, 'xRes': xRes_val, 'dx': dx_ignored},
                              expected_xlim=[xmin_val, xmax_val - expected_dx],
                              expected_xmin=xmin_val, expected_xmax=xmax_val,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_only_xlim(self):
        xlim_val = [-5.0, 21.0]
        expected_xRes = DEFAULT_XRES
        expected_dx = (xlim_val[1] - xlim_val[0]) / (expected_xRes - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx
        self.assertAxisConfig(dim=2, name='x', args={'xlim': xlim_val},
                              expected_xlim=xlim_val,
                              expected_xmin=expected_xmin, expected_xmax=expected_xmax,
                              expected_xRes=expected_xRes, expected_dx=expected_dx)

    def test_xlim_and_dx(self):
        xlim_val = [-5.0, 21.0]
        dx_val = 0.1
        expected_xRes = round((xlim_val[1] - xlim_val[0]) / dx_val) + 1
        expected_dx = (xlim_val[1] - xlim_val[0]) / (expected_xRes - 1) # dx recalculated
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx
        self.assertAxisConfig(dim=2, name='x', args={'xlim': xlim_val, 'dx': dx_val},
                              expected_xlim=xlim_val,
                              expected_xmin=expected_xmin, expected_xmax=expected_xmax,
                              expected_xRes=expected_xRes, expected_dx=expected_dx)

    def test_xlim_and_xRes(self):
        xlim_val = [-5.0, 21.0]
        xRes_val = 56
        expected_dx = (xlim_val[1] - xlim_val[0]) / (xRes_val - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx
        self.assertAxisConfig(dim=2, name='x', args={'xlim': xlim_val, 'xRes': xRes_val},
                              expected_xlim=xlim_val,
                              expected_xmin=expected_xmin, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_xlim_xRes_dx_hierarchy(self):
        # dx should be ignored
        xlim_val = [-5.0, 21.0]
        xRes_val = 56
        dx_ignored = 0.1
        expected_dx = (xlim_val[1] - xlim_val[0]) / (xRes_val - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx
        self.assertAxisConfig(dim=2, name='x', args={'xlim': xlim_val, 'xRes': xRes_val, 'dx': dx_ignored},
                              expected_xlim=xlim_val,
                              expected_xmin=expected_xmin, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    def test_all_hierarchy(self):
        # xlim should take precedence over all others
        xlim_val = [-5.0, 21.0]
        xmin_ignored = -10.0
        xmax_ignored = 30.0
        xRes_val = 56 # This one is still used with xlim
        dx_ignored = 0.1
        expected_dx = (xlim_val[1] - xlim_val[0]) / (xRes_val - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx
        args = {'xlim': xlim_val, 'xmin': xmin_ignored, 'xmax': xmax_ignored, 'xRes': xRes_val, 'dx': dx_ignored}
        self.assertAxisConfig(dim=2, name='x', args=args,
                              expected_xlim=xlim_val,
                              expected_xmin=expected_xmin, expected_xmax=expected_xmax,
                              expected_xRes=xRes_val, expected_dx=expected_dx)

    # --- Tests for dimension-dependent defaults ---

    def test_1D_defaults(self):
        # x-axis (should be standard defaults)
        self.assertAxisConfig(dim=1, name='x', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)
        # y-axis (should have xmax=1, xRes=1)
        self.assertAxisConfig(dim=1, name='y', args={},
                              expected_xlim=[DEFAULT_XMIN, 1.0 - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=1.0,
                              expected_xRes=1, expected_dx=DEFAULT_DX)
        # z-axis (should have xmax=1, xRes=1)
        self.assertAxisConfig(dim=1, name='z', args={},
                              expected_xlim=[DEFAULT_XMIN, 1.0 - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=1.0,
                              expected_xRes=1, expected_dx=DEFAULT_DX)

    def test_2D_defaults(self):
        # x-axis (should be standard defaults)
        self.assertAxisConfig(dim=2, name='x', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)
        # y-axis (should be standard defaults)
        self.assertAxisConfig(dim=2, name='y', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)
        # z-axis (should have xmax=1, xRes=1)
        self.assertAxisConfig(dim=2, name='z', args={},
                              expected_xlim=[DEFAULT_XMIN, 1.0 - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=1.0,
                              expected_xRes=1, expected_dx=DEFAULT_DX)

    def test_3D_defaults(self):
        # x-axis (should be standard defaults)
        self.assertAxisConfig(dim=3, name='x', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)
        # y-axis (should be standard defaults)
        self.assertAxisConfig(dim=3, name='y', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)
        # z-axis (should be standard defaults)
        self.assertAxisConfig(dim=3, name='z', args={},
                              expected_xlim=[DEFAULT_XMIN, DEFAULT_XMAX - DEFAULT_DX],
                              expected_xmin=DEFAULT_XMIN, expected_xmax=DEFAULT_XMAX,
                              expected_xRes=DEFAULT_XRES, expected_dx=DEFAULT_DX)

    # --- Tests using BaseSystem Initialization ---
    def test_base_system_init_only_xlim(self):
        xlim_val = [-5.0, 21.0]
        bs = cf.BaseSystem(1, xlim=xlim_val)
        expected_xRes = DEFAULT_XRES
        expected_dx = (xlim_val[1] - xlim_val[0]) / (expected_xRes - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx

        np.testing.assert_allclose(bs.xlim, xlim_val, rtol=1e-7, atol=1e-7)
        self.assertAlmostEqual(bs.xmin, expected_xmin, places=7)
        self.assertAlmostEqual(bs.xmax, expected_xmax, places=7)
        self.assertEqual(bs.xRes, expected_xRes)
        self.assertAlmostEqual(bs.dx, expected_dx, places=7)

    def test_base_system_init_xmin_xmax_xRes(self):
        xmin_val = -5.0
        xmax_val = 21.0
        xRes_val = 56
        bs = cf.BaseSystem(1, xmin=xmin_val, xmax=xmax_val, xRes=xRes_val)

        expected_dx = (xmax_val - xmin_val) / xRes_val
        expected_xlim = [xmin_val, xmax_val - expected_dx]

        np.testing.assert_allclose(bs.xlim, expected_xlim, rtol=1e-7, atol=1e-7)
        self.assertAlmostEqual(bs.xmin, xmin_val, places=7)
        self.assertAlmostEqual(bs.xmax, xmax_val, places=7)
        self.assertEqual(bs.xRes, xRes_val)
        self.assertAlmostEqual(bs.dx, expected_dx, places=7)

    def test_base_system_init_all_params_xlim_hierarchy(self):
        xlim_val = [-5.0, 21.0]
        xmin_ignored = -10.0
        xmax_ignored = 30.0
        xRes_val = 56
        dx_ignored = 0.1
        bs = cf.BaseSystem(1, xlim=xlim_val, xmin=xmin_ignored, xmax=xmax_ignored, xRes=xRes_val, dx=dx_ignored)

        # Expected values based on xlim and xRes taking precedence
        expected_dx = (xlim_val[1] - xlim_val[0]) / (xRes_val - 1)
        expected_xmin = xlim_val[0]
        expected_xmax = xlim_val[1] + expected_dx

        np.testing.assert_allclose(bs.xlim, xlim_val, rtol=1e-7, atol=1e-7)
        self.assertAlmostEqual(bs.xmin, expected_xmin, places=7)
        self.assertAlmostEqual(bs.xmax, expected_xmax, places=7)
        self.assertEqual(bs.xRes, xRes_val)
        self.assertAlmostEqual(bs.dx, expected_dx, places=7)


if __name__ == '__main__':
    unittest.main()