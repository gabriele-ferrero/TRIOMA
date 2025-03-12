import unittest
import sys
import os
import matplotlib

from TRIOMA.tools.Extractors.PipeSubclasses import CustomTurbulator

matplotlib.use("Agg")
sys.path.append(os.path.abspath("."))
import unittest
from src.TRIOMA.tools.component_tools import (
    Component,
    Fluid,
    Membrane,
    FluidMaterial,
    BreedingBlanket,
    SolidMaterial,
    Geometry,
    Circuit,
    WireCoil,
    GLC,
    GLC_Gas,
)
from src.TRIOMA.tools.materials import Flibe, Steel, Sodium, LiPb
from io import StringIO
from unittest.mock import patch
import pytest
import numpy as np
import matplotlib.pyplot as plt


class TestMSComponent(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some initial values
        self.T_fluid = 300
        self.D_fluid = 1e-9
        self.Solubility_fluid = 0.5
        self.MS_fluid = True
        self.mu_fluid = 1e-3
        self.rho_fluid = 1000
        self.k_fluid = 0.5
        self.cp_fluid = 1.0
        self.k_t_fluid = 0.1
        self.U0_fluid = 0.2
        self.d_Hyd_fluid = 0.3
        self.L_geom = 1.0
        self.thick_geom = 0.5
        self.D_geom = 0.3
        self.k_d_membrane = 1e7
        self.D_membrane = 0.4
        self.thick_membrane = 0.5
        self.K_S_membrane = 0.6
        self.T_membrane = 300
        self.k_r_membrane = 1e7
        self.k_membrane = 0.8
        self.c_in_component = 0.5
        self.eff_component = 0.8
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
        geometry = Geometry(L=1.0, thick=0.5, D=0.3)
        membrane = Membrane(D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8, k_d=1e7)
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.outlet_c_comp()
        self.assertAlmostEqual(
            self.component.c_out, self.c_in_component * (1 - self.eff_component)
        )

    def test_T_leak(self):
        # Test the T_leak() method
        leak = self.component.T_leak()
        self.assertAlmostEqual(
            leak,
            self.c_in_component
            * (self.eff_component)
            * self.D_geom**2
            / 4
            * np.pi
            * self.U0_fluid,
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertEqual(
            self.component.H,
            self.component.membrane.k_d
            / (self.component.fluid.k_t * self.component.fluid.Solubility),
        )
        self.assertEqual(
            self.component.W,
            2
            * self.component.membrane.k_d
            * self.component.geometry.thick
            * (self.component.c_in / self.component.fluid.Solubility) ** 0.5
            / (self.component.membrane.D * self.component.membrane.K_S),
        )

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method ##TODO: Hard coded test
        self.component.use_analytical_efficiency()
        self.assertAlmostEqual(self.component.eff, 0.99871670123992)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency()

        self.assertAlmostEqual(self.component.eff, 0.998984924629)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.analytical_efficiency()
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
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    def test_component_inventory(self):
        self.component.use_analytical_efficiency()
        self.component.get_inventory()
        self.assertEqual(self.component.inv, 0.00531677445914132)

    # def test_inspect(self):
    #     result = "c_in: 0.5\neff: 0.8\nL: 1.0\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 300\n    Solubility: 0.5\n    MS: True\n    D: 1e-09\n    k_t: 0.1\n    d_Hyd: 0.3\n    mu: 0.001\n    rho: 1000\n    U0: 0.2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 300\n    D: 0.4\n    thick: 0.5\n    k_d: 10000000.0\n    K_S: 0.6\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())
    #     result = ""
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect("kjhgfd")
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())


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
        geometry = Geometry(L=1.0, thick=0.5, D=0.3)
        membrane = Membrane(k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8)
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.outlet_c_comp()
        self.assertAlmostEqual(
            self.component.c_out, self.component.c_in * (1 - self.component.eff)
        )

    def test_T_leak(self):
        # Test the T_leak() method
        leak = self.component.T_leak()
        self.assertAlmostEqual(
            leak,
            self.component.c_in
            * self.component.eff
            * self.component.geometry.D**2
            / 4
            * np.pi
            * self.component.fluid.U0,
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(regime, "Mixed regime")

    def test_get_adimensionals(self):
        # Test the get_adimensionals() method
        self.component.get_adimensionals()
        self.assertEqual(
            self.component.H,
            self.component.membrane.k_r
            / self.component.membrane.D
            * self.component.membrane.K_S
            * self.component.membrane.thick
            * self.component.c_in
            / self.component.fluid.Solubility
            * self.component.membrane.D
            * self.component.membrane.K_S
            / (
                self.component.fluid.k_t
                * self.component.fluid.Solubility
                * self.component.membrane.thick
            ),
        )
        self.assertEqual(
            self.component.W,
            self.component.membrane.k_r
            / self.component.membrane.D
            * self.component.membrane.K_S
            * self.component.membrane.thick
            * self.component.c_in
            / self.component.fluid.Solubility,
        )

    def test_use_analytical_efficiency(self):
        # Test the use_analytical_efficiency() method
        self.component.use_analytical_efficiency()
        self.assertAlmostEqual(self.component.eff, 0.998295638580)

    def test_get_efficiency(self):
        # Test the get_efficiency() method
        self.component.get_efficiency(c_guess=self.component.c_in / 2, plotvar=True)
        self.assertAlmostEqual(self.component.eff, 0.9986246, places=5)

    def test_analytical_efficiency(self):
        # Test the analytical_efficiency() method
        self.component.fluid.update_attribute("k_t", None)
        self.component.analytical_efficiency()
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
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )

    # def test_inspect(self):
    #     result = "c_in: 0.5\neff: 0.8\nL: 1.0\nfluid is a <class 'tools.component_tools.Fluid'> class, printing its variables:\n    T: 300\n    Solubility: 0.5\n    MS: False\n    D: 1e-09\n    k_t: 0.1\n    d_Hyd: 0.3\n    mu: 0.001\n    rho: 1000\n    U0: 0.2\n    k: 0.5\n    cp: 1.0\nmembrane is a <class 'tools.component_tools.Membrane'> class, printing its variables:\n    T: 300\n    D: 0.4\n    thick: 0.5\n    k_d: 10000000.0\n    K_S: 0.6\n    k_r: 10000000.0\n    k: 0.8\nH: None\nW: None"
    #     expected_output = result
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

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

    def test_LiPb(self):
        T = 800
        fluid_material = LiPb(T)
        self.component.fluid.set_properties_from_fluid_material(fluid_material)
        self.assertEqual(self.component.fluid.rho, fluid_material.rho)
        solid_material = Steel(T)
        self.component.membrane.set_properties_from_solid_material(solid_material)
        self.assertEqual(self.component.membrane.K_S, solid_material.K_S)

    def test_sodiium(self):
        T = 800
        fluid_material = Sodium(T)
        self.component.fluid.set_properties_from_fluid_material(fluid_material)
        self.assertEqual(self.component.fluid.rho, fluid_material.rho)
        solid_material = Steel(T)
        self.component.membrane.set_properties_from_solid_material(solid_material)
        self.assertEqual(self.component.membrane.K_S, solid_material.K_S)

    def test_update_T_properties(self):
        T = 999
        self.component.update_attribute("T", T)
        self.assertEqual(self.component.fluid.T, T)

    def test_component_inventory(self):
        self.component.get_inventory()
        self.assertEqual(self.component.inv, 0.007008033266553368)
        inv_m_an = self.component.analytical_solid_inventory()
        inv_f_an = self.component.analytical_fluid_inventory()
        inv_m_num = self.component.get_solid_inventory(flag_an=False)
        inv_f_num = self.component.get_fluid_inventory(flag_an=False)
        self.assertAlmostEqual(inv_m_an, inv_m_num)
        self.assertAlmostEqual(inv_f_an, inv_f_num)

    def test_run_inspect(self):
        self.component.inspect()
        self.component.geometry.inspect()

    def test_pumping_power(self):
        self.component.get_pumping_power()
        self.assertAlmostEqual(self.component.pumping_power, 0.019029194163191265)

    def test_plot_component(self):
        self.component.use_analytical_efficiency()
        self.component.outlet_c_comp()
        self.component.plot_component()
        plt.close()

    def test_component_inventory(self):
        self.component.get_inventory(flag_an=False)
        self.assertAlmostEqual(self.component.inv, 0.0070080332665533665, places=7)


class TestExoticComponent(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some initial values
        self.T_fluid = 300
        self.D_fluid = 1e-9
        self.Solubility_fluid = 0.5
        self.MS_fluid = True
        self.mu_fluid = 1e-3
        self.rho_fluid = 1000
        self.k_fluid = 0.5
        self.cp_fluid = 1.0
        self.k_t_fluid = 0.1
        self.U0_fluid = 0.2
        self.d_Hyd_fluid = 0.3
        self.L_geom = 1.0
        self.thick_geom = 0.5
        self.D_geom = 0.3
        self.k_d_membrane = 1e7
        self.D_membrane = 0.4
        self.thick_membrane = 0.5
        self.K_S_membrane = 0.6
        self.T_membrane = 300
        self.k_r_membrane = 1e7
        self.k_membrane = 0.8
        self.c_in_component = 0.5
        self.eff_component = 0.8
        fluid = Fluid(
            T=300,
            D_0=1e-9,
            E_d=0.5,
            Solubility_0=0.5,
            E_s=0.5,
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            k_t=0.1,
            U0=0.2,
            d_Hyd=0.3,
            recirculation=-0.5,
        )
        wirecoil = WireCoil(pitch=1e-2)
        geometry = Geometry(L=1.0, thick=0.5e-3, D=0.3, turbulator=wirecoil)
        membrane = Membrane(
            D_0=0.4,
            E_d=0.1,
            E_S=0.1,
            thick=0.5,
            K_S_0=0.6,
            T=300,
            k_r=1e7,
            k=0.8,
            k_d=1e7,
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_updateTproperties(self):
        self.assertEqual(self.component.fluid.D, 3.984462016634033e-18)
        self.assertEqual(self.component.membrane.D, 0.008358607447235438)
        T = 999
        self.component.update_attribute("T", T)
        self.assertEqual(self.component.fluid.T, T)
        self.assertAlmostEqual(self.component.membrane.D, 0.12519232082575535, places=7)
        self.assertEqual(self.component.fluid.D, 3.003229324075595e-12)

    def test_update_recirculation(self):
        self.assertAlmostEqual(self.component.fluid.recirculation, -0.5, places=5)
        self.assertAlmostEqual(self.component.fluid.U0, 0.1, places=5)

    def test_wirecoil_kt(self):
        self.component.fluid.update_attribute("k_t", None)
        self.component.use_analytical_efficiency()
        self.assertAlmostEqual(self.component.fluid.k_t, 3.74123654043251e-11)
        newturb = CustomTurbulator(a=1, b=1, c=1)
        self.component.geometry.turbulator = newturb
        self.component.fluid.update_attribute("k_t", None)
        self.component.use_analytical_efficiency()
        self.assertAlmostEqual(self.component.fluid.k_t, 0.1)

    def test_cout(self):
        self.component.outlet_c_comp()
        self.assertAlmostEqual(self.component.c_out, 0.3)
        self.component.fluid.update_attribute("recirculation", 0.5)
        self.component.outlet_c_comp()
        self.assertAlmostEqual(self.component.c_out, 5.399292025261743e-07)

    def test_flowrate(self):
        self.component.get_total_flowrate()
        self.assertAlmostEqual(self.component.flowrate, 0.007068583470577035)

    def test_volume(self):
        self.component.define_component_volumes()
        self.assertAlmostEqual(self.component.V, 0.07115628820564542)
        self.assertAlmostEqual(self.component.membrane.V, 0.0004704534998750733)
        self.assertAlmostEqual(self.component.fluid.V, 0.07068583470577035)
        self.component.geometry.get_total_volume()


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

    # def test_inspect(self):
    #     result = "T: 300\nD: 1e-09\nSolubility: 0.5\nMS: True\nmu: 0.001\nrho: 1000\nk: 0.5\ncp: 1.0"
    #     expected_output = result
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.fluid_material.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

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

    # def test_inspect(self):
    #     result = "T: 300\nD: 1e-07\nK_S: 1\nk: None"
    #     expected_output = result
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.solid_material.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), expected_output.strip())

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

    def plot_test(self):
        fig = self.component.plot_component()

    def test_outlet_c_comp(self):
        # Test the outlet_c_comp() method
        self.component.get_cout(print_var=True)
        self.assertAlmostEqual(self.component.c_out, 0.0001473990666223908)

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

    def test_plot_component(self):
        self.component.get_cout()
        self.component.plot_component()
        plt.close()

    # def test_inspect(self):
    #     result = "\n        c_in: 0\nQ: 500000000.0\nTBR: 1.05\nT_out: 900\nT_in: 800\nfluid is a <class 'tools.component_tools.FluidMaterial'> class, printing its variables:\n    T: 850\n    D: 2.4399672021371085e-09\n    Solubility: 0.000454\n    MS: True\n    mu: 0.009616515365587481\n    rho: 1998.2\n    k: 1.1\n    cp: 2386"
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())
    #     result = ""
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect("kjhgfd")
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())


class TestMembrane(unittest.TestCase):
    def setUp(self):
        self.component = Membrane(
            k_d=1e7, D=0.4, thick=0.5, K_S=0.6, T=300, k_r=1e7, k=0.8
        )

    # def test_inspect(self):
    #     result = "T: 300\nD: 0.4\nthick: 0.5\nk_d: 10000000.0\nK_S: 0.6\nk_r: 10000000.0\nk: 0.8"
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())
    #     result = ""
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect("kjhgfd")
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())

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

    # def test_inspect(self):
    #     result = "T: 300\nSolubility: 0.5\nMS: True\nD: 1e-09\nk_t: 0.1\nd_Hyd: 0.3\nmu: 0.001\nrho: 1000\nU0: 0.2\nk: 0.5\ncp: 1.0"
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect()
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())
    #     result = ""
    #     with patch("sys.stdout", new=StringIO()) as fake_out:
    #         self.component.inspect("kjhgfd")
    #         self.assertEqual(fake_out.getvalue().strip(), result.strip())

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
        # result = "k_t is already defined"
        # with patch("sys.stdout", new=StringIO()) as fake_out:
        #     self.component.get_kt()
        #     self.assertEqual(fake_out.getvalue().strip(), result.strip())
        # test when Reynolds number is too low
        self.component.update_attribute("U0", 1e-6)
        self.component.update_attribute("k_t", None)

        self.component.get_kt()
        self.assertEqual(
            (3.66 / self.component.d_Hyd * self.component.D), self.component.k_t
        )
        self.component.update_attribute("U0", 0.008)
        # test when Reynolds number is in a different correlation range
        self.component.get_kt()
        self.assertAlmostEqual(self.component.k_t, 1.2200000000000002e-08)


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
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e7, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Diffusion Limited")

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestMSComponentMixedDiffusionMassTransport(unittest.TestCase):
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
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e7, D=1e-6, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestMSComponentMixedDiffusionSurface(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-2,  # Changed diffusion coefficient
            Solubility=1e-4,  # Changed solubility
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=1,
            d_Hyd=2e-2,
        )
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e-10, D=1e-7, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e-10, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in)


class TestMSComponentMixedMassTransportSurface(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-8,  # Changed diffusion coefficient
            Solubility=1e-4,  # Changed solubility
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=1,
            d_Hyd=2e-2,
        )
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e-8, D=1e-2, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e-8, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


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
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e7, D=1e-2, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mass transport limited")

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3, c_guess=0.3)
        self.assertAlmostEqual(self.component.J_perm, -9.6149466095734e-05, places=5)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestMSComponentSurfaceLimited(unittest.TestCase):
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
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e-16, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e-16, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Surface limited")

    def test_get_efficiency(self):
        # Test that get efficiency works
        self.component.get_efficiency(c_guess=self.component.c_in / 2)


class TestMSComponentFullyMixed(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-11,  # Changed diffusion coefficient
            Solubility=1e-4,  # Changed solubility
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=1,
            d_Hyd=2e-2,
        )
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e-10, D=1e-7, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e-10, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime()
        self.assertEqual(regime, "Mixed regime")

    def test_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


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
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e7, D=1e-2, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Mass transport limited",
        )

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3, c_guess=0.3)
        self.assertAlmostEqual(self.component.J_perm, -4.80747330478670e-05, places=5)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestLMComponentMixedDiffusionMassTransport(unittest.TestCase):
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
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e7, D=1e-7, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Mixed regime",
        )

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestLMComponentMixedDiffusionSurface(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-3,  # Changed diffusion coefficient
            Solubility=1e-2,  # Changed solubility
            MS=False,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=2,
            d_Hyd=2e-3,
        )
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e-3, D=1e-7, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e-3, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Mixed regime",
        )

    def test_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


class TestLMComponentMixedMassTransferSurface(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-6,  # Changed diffusion coefficient
            Solubility=1e-2,  # Changed solubility
            MS=False,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=2,
            d_Hyd=2e-3,
        )
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e-3, D=1e-2, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e-3, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Transport and surface limited regime",
        )

    def test_get_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


class TestLMComponentFullyMixed(unittest.TestCase):
    def setUp(self):
        # Create a Component object with some different initial values
        fluid = Fluid(
            T=750,  # Changed temperature
            D=2e-11,  # Changed diffusion coefficient
            Solubility=1e-2,  # Changed solubility
            MS=False,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            U0=2,
            d_Hyd=2e-3,
        )
        geometry = Geometry(L=1.0, thick=1e-4, D=2e-3)
        membrane = Membrane(
            k_d=1e-9, D=1e-9, thick=1e-4, K_S=0.6e-2, T=700, k_r=1e-9, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(
            regime,
            "Mixed regime",
        )

    def test_efficiency(self):
        # Test the efficiency_vs_analytical() method

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


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
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e7, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e7, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def plot_test(self):
        fig = self.component.plot_component()

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(regime, "Diffusion Limited")

    def test_get_flux(self):
        # Test the get_flux() method
        flux = self.component.get_flux(c=0.3)
        self.assertAlmostEqual(flux, 0.29990976788036605, places=5)

    def test_efficiency_vs_analytical(self):
        # Test the efficiency_vs_analytical() method
        self.component.analytical_efficiency()
        self.component.get_efficiency(c_guess=self.component.c_in / 2)
        self.assertAlmostEqual(
            abs(self.component.eff - self.component.eff_an) / self.component.eff_an,
            0,
            places=2,
        )


class TestLMComponentSurfaceLimited(unittest.TestCase):
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
        geometry = Geometry(L=1.0, thick=1e-2, D=2e-2)
        membrane = Membrane(
            k_d=1e-16, D=1e-9, thick=1e-2, K_S=0.6e-2, T=700, k_r=1e-16, k=0.8
        )
        self.component = Component(
            c_in=0.5, geometry=geometry, eff=0.8, fluid=fluid, membrane=membrane
        )

    def plot_test(self):
        self.component.plot()

    def test_get_regime(self):
        regime = self.component.get_regime(print_var=True)
        self.assertEqual(regime, "Surface limited")

    def test_get_efficiency(self):
        # Test that get efficiency works

        self.component.get_efficiency(c_guess=self.component.c_in / 2)


class testclosedCircuit(unittest.TestCase):
    def setUp(self):
        fluid = Fluid(
            T=900,
            D=1e-7,
            Solubility=0.5,
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            k_t=0.1,
            U0=1,
            d_Hyd=2e-2,
        )
        geometry = Geometry(L=10.0, thick=0.5e-3, D=2e-2)
        membrane = Membrane(
            D_0=1e-7, E_d=1, thick=0.5e-3, K_S=0.6, T=900, k_r=1e7, k=0.8, k_d=1e7
        )
        component = Component(
            geometry=geometry, fluid=fluid, membrane=membrane, loss=False
        )
        component2 = Component(
            geometry=geometry, fluid=fluid, membrane=membrane, loss=False, name="PAV"
        )
        componentBB = BreedingBlanket(
            c_in=1e-3,
            Q=0.5e9,
            TBR=1.05,
            T_out=900,
            T_in=800,
            fluid=Flibe(850),
            name="BB",
        )
        geometry_HX = Geometry(L=10.0, thick=1e-3, D=2e-3)
        membrane_HX = Membrane(
            D_0=1e-9, E_d=0.2, thick=1e-3, K_S=0.6, T=850, k_r=1e7, k=0.8, k_d=1e7
        )
        fluid_HX = Fluid(
            T=850,
            D=1e-9,
            Solubility=0.5,
            MS=True,
            mu=1e-3,
            rho=1000,
            k=0.5,
            cp=1.0,
            k_t=0.1,
            U0=1,
            d_Hyd=2e-2,
        )
        self.component_HX = Component(
            c_in=1e-3,
            geometry=geometry_HX,
            fluid=fluid_HX,
            membrane=membrane_HX,
            loss=True,
        )
        HX = self.component_HX.split_HX(
            N=11,
            T_in_hot=900,
            T_out_hot=800,
            T_in_cold=581,
            T_out_cold=800,
            R_sec=0,
            Q=0.5e9 / 5e3,
            plotvar=True,
        )
        plt.close()
        self.circuit = Circuit([componentBB, component, HX, component2], closed=True)

    def test_circuit(self):
        # Only test that Circuit functions are called and return correctly
        # Test closed circuit
        self.circuit.solve_circuit()
        self.circuit.get_eff_circuit()
        self.circuit.get_gains_and_losses()
        self.circuit.inspect_circuit()
        # fig=self.circuit.plot_circuit()
        # Test component split and then test open circuit
        self.component_HX.converge_split_HX(
            T_in_hot=900,
            T_out_hot=800,
            T_in_cold=581,
            T_out_cold=800,
            R_sec=0,
            Q=0.5e9 / 5e3,
            plotvar=True,
        )
        Circuit_HX = self.component_HX.split_HX(
            N=11,
            T_in_hot=900,
            T_out_hot=800,
            T_in_cold=581,
            T_out_cold=800,
            R_sec=0,
            Q=0.5e9 / 5e3,
        )
        Circuit_HX.solve_circuit()
        Circuit_HX.inspect_circuit()
        Circuit_HX.add_component(self.component_HX)
        Circuit_HX.solve_circuit()
        Circuit_HX.plot_circuit()
        metal_cost = np.ones(len(Circuit_HX.components))
        fluid_cost = np.ones(len(Circuit_HX.components))
        Circuit_HX.estimate_cost(metal_cost=metal_cost, fluid_cost=fluid_cost)
        plt.close()
        self.circuit.plot_circuit()

        plt.close()
        Circuit_HX.get_inventory()
        self.circuit.inspect_circuit(name="PAV")
        self.circuit.get_circuit_pumping_power()


class testLMGLCComponent(unittest.TestCase):
    def setUp(self):
        R_const = 8.314
        T = 400 + 273.15
        c_in = 1e-2
        c_out = 9e-3
        Z = 0.6
        R = 0.3
        SweepGas = GLC_Gas(G_gas=3 * 1e-3 / 3600, pg_in=0, p_tot=1.5e5)
        Q_l = 71.0 * 1e-3 / 3600
        Flibe = Fluid(Solubility=1.33e-4 * np.exp(-1350 / R_const / T), MS=False)
        Melodie = GLC(
            H=Z,
            R=R,
            c_in=c_in,
            fluid=Flibe,
            GLC_gas=SweepGas,
            T=T,
            G_L=Q_l,
            c_out=c_out,
        )
        self.GLC = Melodie

    def test_get_kla(self):
        self.GLC.get_kla_from_cout()
        self.assertAlmostEqual(self.GLC.kla, 1.2850594291115214e-05, places=8)

    def test_get_cout(self):
        self.GLC.kla = 1.2850594291115214e-05
        self.GLC.get_c_out()

        self.assertAlmostEqual(self.GLC.c_out, 0.009, places=7)

    def test_get_z_from_eff(self):
        self.GLC.kla = 1.2850594291115214e-05
        self.GLC.get_z_from_eff()
        self.assertAlmostEqual(self.GLC.H, 0.6, places=7)


class testMSGLCComponent(unittest.TestCase):
    def setUp(self):
        R_const = 8.314
        T = 400 + 273.15
        c_in = 1e-2
        c_out = 9e-3
        Z = 0.6
        R = 0.3
        SweepGas = GLC_Gas(G_gas=3 * 1e-3 / 3600, pg_in=0, p_tot=1.5e5)
        Q_l = 71.0 * 1e-3 / 3600
        Flibe = Fluid(Solubility=1.33e-4 * np.exp(-1350 / R_const / T), MS=True)
        Melodie = GLC(
            H=Z,
            R=R,
            c_in=c_in,
            fluid=Flibe,
            GLC_gas=SweepGas,
            T=T,
            G_L=Q_l,
            c_out=c_out,
        )
        self.GLC = Melodie

    def test_get_kla(self):
        self.GLC.get_kla_from_cout()
        self.assertAlmostEqual(self.GLC.kla, 2.9765872207306292e-05, places=8)

    def test_get_cout(self):
        self.GLC.kla = 2.9765872207306292e-05
        self.GLC.get_c_out()

        self.assertAlmostEqual(self.GLC.c_out, 0.009, places=7)

    def test_get_z_from_eff(self):
        self.GLC.kla = 2.9765872207306292e-05
        self.GLC.get_z_from_eff()
        self.assertAlmostEqual(self.GLC.H, 0.6, places=7)


if __name__ == "__main__":
    unittest.main()
