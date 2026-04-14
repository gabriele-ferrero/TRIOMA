"""Component tools and utilities for TRIOMA."""

from typing import Union, Optional, TYPE_CHECKING

from TRIOMA.tools.TriomaClass import TriomaClass
from TRIOMA.tools.Extractors.PipeSubclasses import (
    Geometry,
    Fluid,
    Membrane,
    FluidMaterial,
    SolidMaterial,
    Turbulator,
    WireCoil,
    CustomTurbulator,
)
from TRIOMA.tools.Extractors.GasLiquidContactor import GLC

__all__ = [
    "TriomaClass",
    "Geometry",
    "Fluid",
    "Membrane",
    "FluidMaterial",
    "SolidMaterial",
    "Turbulator",
    "WireCoil",
    "CustomTurbulator",
    "GLC",
]


if TYPE_CHECKING:
    from TRIOMA.tools.Extractors.PAV import Component
    from TRIOMA.tools.BreedingBlanket import BreedingBlanket


def connect_to_component(
    self: TriomaClass, component2: Optional[Union["Component", "BreedingBlanket", "GLC"]] = None
) -> None:
    """
    Connect this component to another component.

    Sets the inlet concentration of component2 equal to the outlet of self.

    Args:
        self: The TRIOMA component instance.
        component2: The component to connect to.
    """
    if component2 is None:
        raise ValueError("component2 cannot be None")
    component2.update_attribute("c_in", self.c_out)


# Dynamically add connect_to_component method to component classes
from TRIOMA.tools.Extractors.PAV import Component
from TRIOMA.tools.BreedingBlanket import BreedingBlanket

Component.connect_to_component = connect_to_component
BreedingBlanket.connect_to_component = connect_to_component
GLC.connect_to_component = connect_to_component
