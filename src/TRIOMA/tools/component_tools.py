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
from TRIOMA.tools.Extractors.GasLiquidContactor import GLC_Gas, GLC

import TRIOMA.tools.correlations as corr
import matplotlib.pyplot as plt

from typing import Union
import matplotlib.pyplot as plt
from scipy import integrate
from TRIOMA.tools.Extractors import extractor
from TRIOMA.tools.Extractors.PAV import Component
from TRIOMA.tools.BreedingBlanket import BreedingBlanket
from TRIOMA.tools.Circuit import Circuit
import numpy as np


def connect_to_component(
    self, component2: Union["Component", "BreedingBlanket", "GLC"] = None
):
    """sets the inlet conc of the object component equal to the outlet of self"""
    component2.update_attribute("c_in", self.c_out)


# Component.converge_split_HX = converge_split_HX
# Component.split_HX = split_HX
Component.connect_to_component = connect_to_component
BreedingBlanket.connect_to_component = connect_to_component
GLC.connect_to_component = connect_to_component
