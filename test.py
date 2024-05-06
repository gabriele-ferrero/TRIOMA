import unittest
from tools.component_tools import FluidMaterial


class FluidMaterialTests(unittest.TestCase):
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
        expected_output = "T: 300\nD: 1e-09\nSolubility: 0.5\nMS: 28.97\nmu: 0.001\nrho: 1000\nk: 0.5\ncp: 1.0\n"
        self.assertEqual(self.fluid_material.inspect(), expected_output)


if __name__ == "__main__":
    unittest.main()
