import unittest
from tools.component_tools import Component, Fluid, Membrane


class TestComponent(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some initial values
        fluid = Fluid(
            T=300,
            D=1e-9,
            Solubility=0.5,
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            k_t=0.1,
            U0=0.2,
            d_Hyd=0.3,
        )
        membrane = Membrane(k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8)
        self.component = Component(c_in=0.5, eff=0.8, fluid=fluid, membrane=membrane)

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.outlet_c_comp()
        self.assertAlmostEqual(self.component.c_out, 0.1)

    def test_T_leak(self):
        # Test the T_leak() method
        leak = self.component.T_leak()
        self.assertAlmostEqual(leak, 0.4)

    def test_get_regime(self):
        self.component2 = Component(c_in=0.5, eff=0.8)
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No fluid selected")
        self.component2.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        # Test the get_regime() method
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No membrane selected")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method

        self.component.get_adimensionals()
        self.assertEqual(self.component.H, 200000000)
        self.assertEqual(self.component.W, 41666666.66666667)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.99871669093034)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)

        self.assertAlmostEqual(self.component.eff, 0.998984924629, places=2)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.9987166909303, places=2)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.0015, places=2)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 2.9215784663, places=4)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


if __name__ == "__main__":
    unittest.main()
