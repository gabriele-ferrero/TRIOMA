import unittest
from tools.component_tools import Component, Fluid, Membrane


class TestComponent(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some initial values
        self.component = Component(c_in=0.5, eff=0.8)

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.outlet_c_comp()
        self.assertEqual(self.component.c_out, 0.1)

    def test_T_leak(self):
        # Test the T_leak() method
        leak = self.component.T_leak()
        self.assertEqual(leak, 0.4)

    def test_get_regime(self):
        # Test the get_regime() method
        self.component.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        self.component.membrane = Membrane(k_d=0.3, D=0.4, thick=0.5, K_S=0.6)
        regime = self.component.get_regime()
        self.assertEqual(regime, "MS regime")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        self.component.membrane = Membrane(k_d=0.3, D=0.4, thick=0.5, K_S=0.6)
        self.component.get_adimensionals()
        self.assertEqual(self.component.H, 0.08)
        self.assertEqual(self.component.W, 0.12)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.analytical_efficiency = lambda L: 0.9
        self.component.use_analytical_efficiency(L=1.0)
        self.assertEqual(self.component.eff, 0.9)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        self.component.J_perm = 0.05
        efficiency = self.component.get_efficiency(L=1.0)
        self.assertAlmostEqual(efficiency, 0.6, places=2)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        self.component.membrane = Membrane(k_d=0.3, D=0.4, thick=0.5, K_S=0.6)
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.7, places=2)

    def test_get_flux(self):
        # Test the get_flux() method
        self.component.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        self.component.membrane = Membrane(k_d=0.3, D=0.4, thick=0.5, K_S=0.6)
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.15, places=2)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        self.component.fluid = Fluid(d_Hyd=0.1, k=0.2, rho=0.3, U0=0.4, mu=0.5, cp=0.6)
        self.component.membrane = Membrane(thick=0.7)
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(U, 0.2857, places=4)


if __name__ == "__main__":
    unittest.main()
