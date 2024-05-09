import unittest
from tools.component_tools import (
    Component,
    Fluid,
    Membrane,
    FluidMaterial,
    BreedingBlanket,
    SolidMaterial,
)
from tools.materials import Flibe, Steel
from io import StringIO
from unittest.mock import patch
import pytest
from tools.materials import Flibe


class TestMSComponent(unittest.TestCase):
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
        regime = self.component2.get_regime(print_var=True)
        self.assertEqual(regime, "No fluid selected")
        self.component2.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        # Test the get_regime() method
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No membrane selected")
        self.component2.membrane = Membrane(
            k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8
        )
        regime = self.component2.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method

        self.component.get_adimensionals()
        self.assertEqual(self.component.H, 200000000)
        self.assertEqual(self.component.W, 41666666.66666667)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.99871670123992)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0)

        self.assertAlmostEqual(self.component.eff, 0.998984924629)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.99871670123992)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3, c_guess=0.3)
        self.assertAlmostEqual(flux, 0.0014967, places=5)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 2.9215784663)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "c_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 300\n    Solubility: 0.5\n    MS: True\n    D: 1e-09\n    k_t: 0.1\n    d_Hyd: 0.3\n    mu: 0.001\n    rho: 1000\n    U0: 0.2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 300\n    D: 0.4\n    thick: 0.5\n    k_d: 10000000.0\n    K_S: 0.6\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestLMComponent(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some initial values
        fluid = Fluid(
            T=300,
            D=1e-9,
            Solubility=0.5,
            MS=False,
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
        regime = self.component2.get_regime(print_var=True)
        self.assertEqual(regime, "No fluid selected")
        self.component2.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=False)
        # Test the get_regime() method
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No membrane selected")
        self.component2.membrane = Membrane(
            k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8
        )
        regime = self.component2.get_regime(print_var=True)
        self.assertEqual(regime, "Mass transport limited")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method

        self.component.get_adimensionals()
        self.assertEqual(self.component.H, 72000000)
        self.assertEqual(self.component.W, 7500000)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.998295638580)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(
            L=1.0, c_guess=self.component.c_in / 2, plotvar=True
        )

        self.assertAlmostEqual(self.component.eff, 0.9986246, places=5)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.fluid.update_attribute("k_t", None)
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.00053628139636452)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.01314458525095)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 2.9215784663)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "c_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 300\n    Solubility: 0.5\n    MS: False\n    D: 1e-09\n    k_t: 0.1\n    d_Hyd: 0.3\n    mu: 0.001\n    rho: 1000\n    U0: 0.2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 300\n    D: 0.4\n    thick: 0.5\n    k_d: 10000000.0\n    K_S: 0.6\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        expected_output = result
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

    def test_set_attribute(self):
        self.component.update_attribute("c_in", 0.6)
        self.assertEqual(self.component.c_in, 0.6)
        ##test internal attribute
        self.component.update_attribute("U0", 0.3)
        self.assertEqual(self.component.fluid.U0, 0.3)
        with pytest.raises(ValueError) as excinfo:
            self.component.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of Component"

    def test_set_material_properties(self):
        T = 800
        fluid_material = Flibe(T)
        self.component.fluid.set_properties_from_fluid_material(fluid_material)
        self.assertEqual(self.component.fluid.rho, fluid_material.rho)
        solid_material = Steel(T)
        self.component.membrane.set_properties_from_solid_material(solid_material)
        self.assertEqual(self.component.membrane.K_S, solid_material.K_S)


class TestFluidMaterial(unittest.TestCase):
    def setUp(self):
        self.fluid_material = FluidMaterial(
            T=300, D=1e-9, Solubility=0.5, MS=True, mu=1e-3, rho=1000, k=0.5, cp=1.0
        )

    def test_attributes(self):
        self.assertEqual(self.fluid_material.T, 300)
        self.assertEqual(self.fluid_material.D, 1e-9)
        self.assertEqual(self.fluid_material.Solubility, 0.5)
        self.assertEqual(self.fluid_material.MS, True)
        self.assertEqual(self.fluid_material.mu, 1e-3)
        self.assertEqual(self.fluid_material.rho, 1000)
        self.assertEqual(self.fluid_material.k, 0.5)
        self.assertEqual(self.fluid_material.cp, 1.0)

    def test_inspect(self):
        result = "T: 300\nD: 1e-09\nSolubility: 0.5\nMS: True\nmu: 0.001\nrho: 1000\nk: 0.5\ncp: 1.0"
        expected_output = result
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.fluid_material.inspect()
            self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

    def test_update_attribute(self):
        self.fluid_material.update_attribute("T", 400)
        self.assertEqual(self.fluid_material.T, 400)
        with pytest.raises(ValueError) as excinfo:
            self.fluid_material.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of FluidMaterial"


class Test_SolidMaterial(unittest.TestCase):
    def setUp(self):
        self.solid_material = SolidMaterial(T=300, K_S=1, D=1e-7)

    def test_attributes(self):
        self.assertEqual(self.solid_material.T, 300)
        self.assertEqual(self.solid_material.K_S, 1)

    def test_inspect(self):
        result = "T: 300\nD: 1e-07\nK_S: 1\nk: None"
        expected_output = result
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.solid_material.inspect()
            self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

    def test_update_attribute(self):
        self.solid_material.update_attribute("T", 400)
        self.assertEqual(self.solid_material.T, 400)
        with pytest.raises(ValueError) as excinfo:
            self.solid_material.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of SolidMaterial"


class Test_BB_Component(unittest.TestCase):
    def setUp(self):
        self.component = BreedingBlanket(
            c_in=0, Q=0.5e9, TBR=1.05, T_out=900, T_in=800, fluid=Flibe(850)
        )

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.get_cout(print_var=True)
        self.assertAlmostEqual(self.component.c_out, 0.0002947981332)

    def test_get_flowrate(self):
        self.component.get_flowrate()
        self.assertAlmostEqual(self.component.m_coolant, 2095.5574182607)

    def test_set_attribute(self):
        self.component.update_attribute("c_in", 0.6)
        self.assertEqual(self.component.c_in, 0.6)
        ##test internal attribute
        self.component.update_attribute("cp", 0.3)
        self.assertEqual(self.component.fluid.cp, 0.3)
        with pytest.raises(ValueError) as excinfo:
            self.component.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of BreedingBlanket"

    def test_inspect(self):
        result = "\n        c_in: 0\nQ: 500000000.0\nTBR: 1.05\nT_out: 900\nT_in: 800\nfluid is a <class 'tools.component_tools.FluidMaterial'> class, printing its variables:\n    T: 850\n    D: 2.4399672021371085e-09\n    Solubility: 0.000454\n    MS: True\n    mu: 0.009616515365587481\n    rho: 1998.2\n    k: 1.1\n    cp: 2386"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestMembrane(unittest.TestCase):
    def setUp(self):
        self.component = Membrane(
            k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8
        )

    def test_inspect(self):
        result = "T: 300\nD: 0.4\nthick: 0.5\nk_d: 10000000.0\nK_S: 0.6\nk_r: 10000000.0\nk: 0.8"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())

    def test_set_attribute(self):
        self.component.update_attribute("thick", 1e-3)
        self.assertEqual(self.component.thick, 1e-3)
        with pytest.raises(ValueError) as excinfo:
            self.component.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of Membrane"


class TestFluid(unittest.TestCase):
    def setUp(self):
        self.component = Fluid(
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

    def test_inspect(self):
        result = "T: 300\nSolubility: 0.5\nMS: True\nD: 1e-09\nk_t: 0.1\nd_Hyd: 0.3\nmu: 0.001\nrho: 1000\nU0: 0.2\nk: 0.5\ncp: 1.0"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())

    def test_set_attribute(self):
        self.component.update_attribute("Solubility", 1e-3)
        self.assertEqual(self.component.Solubility, 1e-3)
        with pytest.raises(ValueError) as excinfo:
            self.component.update_attribute("kghufh", 0.3)
        assert str(excinfo.value) == "'kghufh' is not an attribute of Fluid"

    def test_get_kt(self):
        # Test no d_hyd defined
        self.component.update_attribute("d_Hyd", None)
        result = "Hydraulic Diameter is not defined"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.get_kt()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        # test get_kt
        self.component.update_attribute("d_Hyd", 0.3)
        self.component.update_attribute("k_t", None)
        self.component.get_kt()
        self.assertEqual(self.component.k_t, 8.046408367835323e-06)
        # test when k_t is already defined
        result = "k_t is already defined"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.get_kt()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        # test when Reynolds number is too low
        self.component.update_attribute("U0", 0)
        self.component.update_attribute("k_t", None)
        with pytest.raises(ValueError) as excinfo:
            self.component.get_kt()
        assert str(excinfo.value) == "Reynolds number is too low"
        self.component.update_attribute("U0", 0.008)
        # test when Reynolds number is in a different correlation range
        self.component.get_kt()
        self.assertAlmostEqual(self.component.k_t, 4.25875681945274e-07)


class TestMSComponentDiffusionLimited(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-9,  # Changed diffusion coefficient
            Solubility=1e-4,  # Changed solubility
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=1,
            d_Hyd=2e-2,
        )
        membrane = Membrane(
            k_d=1e7, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
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
        regime = self.component2.get_regime(print_var=True)
        self.assertEqual(regime, "No fluid selected")
        self.component2.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        # Test the get_regime() method
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No membrane selected")
        self.component2.membrane = Membrane(
            k_d=1e10, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e10, k=0.8
        )
        regime = self.component2.get_regime()
        self.assertEqual(regime, "Mixed regime")
        regime = self.component.get_regime()
        self.assertEqual(regime, "Diffusion Limited")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertAlmostEqual(self.component.H, 1435637385922034.5)
        self.assertAlmostEqual(self.component.W, 2.3570226039551585e18)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 1.2241635894638e-05)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0)

        self.assertAlmostEqual(self.component.eff, 1.22362589068947e-05)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 1.2241635894638e-05)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.299659862245049, places=5)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 5.3125160576403)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "\n        c_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 750\n    Solubility: 0.0001\n    MS: True\n    D: 2e-09\n    k_t: None\n    d_Hyd: 0.02\n    mu: 0.001\n    rho: 1000\n    U0: 1\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 700\n    D: 1e-09\n    thick: 0.01\n    k_d: 10000000.0\n    K_S: 0.006\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestMSComponentMassTransportLimited(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-9,  # Changed diffusion coefficient
            Solubility=1e-2,  # Changed solubility
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=2,
            d_Hyd=2e-3,
        )
        membrane = Membrane(
            k_d=1e7, D=1e-2, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
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
        regime = self.component2.get_regime(print_var=True)
        self.assertEqual(regime, "No fluid selected")
        self.component2.fluid = Fluid(k_t=0.1, Solubility=0.2, MS=True)
        # Test the get_regime() method
        regime = self.component2.get_regime()
        self.assertEqual(regime, "No membrane selected")
        self.component2.membrane = Membrane(
            k_d=1e10, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e10, k=0.8
        )
        regime = self.component2.get_regime()
        self.assertEqual(regime, "Mixed regime")
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mass transport limited")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertAlmostEqual(self.component.H, 6240284261200.083, places=1)
        self.assertAlmostEqual(self.component.W, 235702260.39551586)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.148068462270955)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0)

        self.assertAlmostEqual(self.component.eff, 0.1481790664501431)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.148068462270955)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3, c_guess=0.3)
        self.assertAlmostEqual(self.component.J_perm, -9.6149466095734e-05, places=5)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 8.892858063573009)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "\nc_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 750\n    Solubility: 0.01\n    MS: True\n    D: 2e-09\n    k_t: None\n    d_Hyd: 0.002\n    mu: 0.001\n    rho: 1000\n    U0: 2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 700\n    D: 0.01\n    thick: 0.0001\n    k_d: 10000000.0\n    K_S: 0.006\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestLMComponentMassTransportLimited(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-9,  # Changed diffusion coefficient
            Solubility=1e-2,  # Changed solubility
            MS=False,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=2,
            d_Hyd=2e-3,
        )
        membrane = Membrane(
            k_d=1e7, D=1e-2, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
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
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Mass transport limited",
        )

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertAlmostEqual(self.component.H, 11232511670.1601, places=1)
        self.assertAlmostEqual(self.component.W, 30000.0)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.1480681147484217)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0)

        self.assertAlmostEqual(self.component.eff, 0.1481790664501431)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.148068114748421)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3, c_guess=0.3)
        self.assertAlmostEqual(self.component.J_perm, -4.80747330478670e-05, places=5)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 8.892858063573009)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "\nc_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 750\n    Solubility: 0.01\n    MS: False\n    D: 2e-09\n    k_t: None\n    d_Hyd: 0.002\n    mu: 0.001\n    rho: 1000\n    U0: 2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 700\n    D: 0.01\n    thick: 0.0001\n    k_d: 10000000.0\n    K_S: 0.006\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestLMComponentDiffusionLimited(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-5,  # Changed diffusion coefficient
            Solubility=1e-4,  # Changed solubility
            MS=False,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=1,
            d_Hyd=2e-2,
        )
        membrane = Membrane(
            k_d=1e7, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
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
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(regime, "Diffusion Limited")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertAlmostEqual(self.component.H, 625629566141.5754)
        self.assertAlmostEqual(self.component.W, 3000000000000000.0)

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff, 0.0017292165186562)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(L=1.0)

        self.assertAlmostEqual(self.component.eff, 0.001729231620799853)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency(L=1.0)
        self.assertAlmostEqual(self.component.eff_an, 0.0017292165186562)

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.299909767880366, places=5)

    def test_get_global_HX_coeff(self):
        # Test the get_global_HX_coeff() method
        U = self.component.get_global_HX_coeff(R_conv_sec=0.1)
        self.assertAlmostEqual(self.component.U, 5.3125160576403)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency(L=1.0)
        self.component.get_efficiency(L=1.0, c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_inspect(self):
        result = "\n        c_in: 0.5\neff: 0.8\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 750\n    Solubility: 0.0001\n    MS: False\n    D: 2e-05\n    k_t: None\n    d_Hyd: 0.02\n    mu: 0.001\n    rho: 1000\n    U0: 1\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 700\n    D: 1e-09\n    thick: 0.01\n    k_d: 10000000.0\n    K_S: 0.006\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect()
            self.assertEqual(fake_out.getvalue().strip(), result.strip())
        result = ""
        with patch("sys.stdout", new=StringIO()) as fake_out:
            self.component.inspect("kjhgfd")
            self.assertEqual(fake_out.getvalue().strip(), result.strip())


if __name__ == "__main__":
    unittest.main()
