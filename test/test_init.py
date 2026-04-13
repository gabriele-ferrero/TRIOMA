"""
Tests for package __init__.py files.

Verifies that the package structure enables proper imports and functioning.
"""

import pytest


class TestRootInit:
    """Test TRIOMA root package __init__.py"""

    def test_import_from_trioma_root(self):
        """Test that main classes are importable from TRIOMA"""
        from TRIOMA import Circuit, BreedingBlanket

        assert Circuit is not None
        assert BreedingBlanket is not None


class TestToolsInit:
    """Test TRIOMA.tools package __init__.py"""

    def test_import_from_tools(self):
        """Test that classes are importable from TRIOMA.tools"""
        from TRIOMA.tools import TriomaClass, Circuit, BreedingBlanket

        assert TriomaClass is not None
        assert Circuit is not None
        assert BreedingBlanket is not None


class TestExtractorsInit:
    """Test TRIOMA.tools.Extractors package __init__.py"""

    def test_import_from_extractors(self):
        """Test that extractor classes are importable from TRIOMA.tools.Extractors"""
        from TRIOMA.tools.Extractors import (
            Component,
            GLC,
            GLC_Gas,
            Geometry,
            Fluid,
            Membrane,
            FluidMaterial,
            SolidMaterial,
            Turbulator,
            WireCoil,
            CustomTurbulator,
        )

        assert Component is not None
        assert GLC is not None
        assert GLC_Gas is not None
        assert Geometry is not None
        assert Fluid is not None
        assert Membrane is not None
        assert FluidMaterial is not None
        assert SolidMaterial is not None
        assert Turbulator is not None
        assert WireCoil is not None
        assert CustomTurbulator is not None


class TestImportConsistency:
    """Test that imports are consistent across different import paths"""

    def test_circuit_same_from_different_paths(self):
        """Test that Circuit imported from different paths is the same class"""
        from TRIOMA import Circuit as Circuit1
        from TRIOMA.tools import Circuit as Circuit2
        from TRIOMA.tools.Circuit import Circuit as Circuit3

        assert Circuit1 is Circuit2
        assert Circuit2 is Circuit3

    def test_component_same_from_different_paths(self):
        """Test that Component imported from different paths is the same class"""
        from TRIOMA.tools.Extractors import Component as Comp1
        from TRIOMA.tools.Extractors.PAV import Component as Comp2

        assert Comp1 is Comp2

    def test_glc_same_from_different_paths(self):
        """Test that GLC imported from different paths is the same class"""
        from TRIOMA.tools.Extractors import GLC as GLC1
        from TRIOMA.tools.Extractors.GasLiquidContactor import GLC as GLC2

        assert GLC1 is GLC2


class TestBackwardCompatibility:
    """Test that old import paths still work for backward compatibility"""

    def test_import_from_pav_directly(self):
        """Test that direct import from PAV still works"""
        from TRIOMA.tools.Extractors.PAV import Component

        assert Component is not None

    def test_import_from_glc_directly(self):
        """Test that direct import from GasLiquidContactor still works"""
        from TRIOMA.tools.Extractors.GasLiquidContactor import GLC, GLC_Gas

        assert GLC is not None
        assert GLC_Gas is not None

    def test_import_from_pipe_subclasses_directly(self):
        """Test that direct import from PipeSubclasses still works"""
        from TRIOMA.tools.Extractors.PipeSubclasses import (
            Geometry,
            Fluid,
            Membrane,
        )

        assert Geometry is not None
        assert Fluid is not None
        assert Membrane is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
