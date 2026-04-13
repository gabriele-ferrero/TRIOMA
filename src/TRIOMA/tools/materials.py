"""Material property functions for TRIOMA."""

import numpy as np

from TRIOMA.tools.component_tools import FluidMaterial, SolidMaterial

# Physical constants
N_A: float = 6.022e23
k_b: float = 8.617e-5
R_const: float = 8.314


def Flibe(T: float) -> FluidMaterial:
    """
    Return Flibe (LiF-BeF2) fluid material properties.

    Args:
        T: Temperature in Kelvin

    Returns:
        FluidMaterial object with Flibe properties
    """

    def density(T: float) -> float:
        return 2413 - 0.488 * T

    def viscosity(T: float) -> float:
        return 1.16e-4 * np.exp(3755 / T)

    def H_diff(T: float) -> float:
        R_const_local = 8.314
        return 9.3e-7 * np.exp(-42e3 / (R_const_local * T))

    def k() -> float:
        # thermal conductivity W/m/K
        return 1.1

    def cp() -> float:
        # specific heat J/kg/K
        return 2386

    def k_H(T: float) -> float:
        # Henry's constant mol/m^3/Pa
        return 4.54e-4

    flibe_material = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_H(T),
        MS=True,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )

    return flibe_material


def Sodium(T: float) -> FluidMaterial:
    """
    Return Sodium liquid metal material properties.

    Args:
        T: Temperature in Kelvin

    Returns:
        FluidMaterial object with Sodium properties
    """

    def density(T: float) -> float:
        return 219 + 275.32 * (1 - T / 2504.7) + 511.58 * (1 - T / 2503.7) ** 0.5

    def viscosity(T: float) -> float:
        return np.exp(-6.4406 - 0.3958 * np.log(T) + 556.835 / T)

    def H_diff(T: float) -> float:
        R_const_local = 8.314
        return 2e-5 * np.exp(-49053 / (R_const_local * T))

    def k() -> float:
        # thermal conductivity W/m/K
        return 1  # TODO

    def cp() -> float:
        # specific heat J/kg/K
        return 1  # TODO

    def k_S(T: float) -> float:
        # Henry's constant mol/m^3/Pa
        return 10 ** (0.86 - 122 / T)

    sodium_material = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_S(T),
        MS=False,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )
    return sodium_material


def LiPb(T: float) -> FluidMaterial:
    """
    Return LiPb (Lithium-Lead) liquid metal material properties.

    Args:
        T: Temperature in Kelvin

    Returns:
        FluidMaterial object with LiPb properties
    """

    def density(T: float) -> float:
        return 9659.8  # TODO

    def viscosity(T: float) -> float:
        return 1  # TODO

    def H_diff(T: float) -> float:
        return 1  # TODO

    def k() -> float:
        # thermal conductivity W/m/K
        return 1  # TODO

    def cp() -> float:
        # specific heat J/kg/K
        return 1  # TODO

    def k_S(T: float) -> float:
        # Henry's constant mol/m^3/Pa
        MM_LiPb = 180  # TODO
        return 4.7e-7 * np.exp(-9e3 / (R_const * T) * density(T) / MM_LiPb)  # TODO

    lipb_material = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_S(T),
        MS=False,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )
    return lipb_material


def Steel(T: float) -> SolidMaterial:
    """
    Return Steel solid material properties.

    Args:
        T: Temperature in Kelvin

    Returns:
        SolidMaterial object with Steel properties
    """

    def H_diff(T: float) -> float:
        D_met = 5.81e-7 * np.exp(-66.3e3 / (R_const * T))
        return D_met

    def K_S(T: float) -> float:
        return 1

    steel_material = SolidMaterial(T=T, D=H_diff(T), K_S=K_S(T))

    return steel_material

