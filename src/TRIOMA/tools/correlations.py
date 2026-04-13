"""Correlations for heat and mass transfer calculations."""

import numpy as np


def Nu_SiederTate(Re: float, Pr: float, mu_w: float, mu_c: float) -> float:
    """
    Calculate the Nusselt number using the Sieder-Tate correlation.

    Args:
        Re: Reynolds number
        Pr: Prandtl number
        mu_w: Viscosity of the wall fluid
        mu_c: Viscosity of the core fluid

    Returns:
        Nusselt number
    """
    return 0.027 * Re ** (4 / 5) * Pr ** (0.3) * (mu_c / mu_w) ** (0.14)


def Nu_Gnielinsky(Re: float, Pr: float, f: float) -> float:
    """
    Calculate the Nusselt number using the Gnielinsky correlation.

    Args:
        Re: Reynolds number
        Pr: Prandtl number
        f: Friction factor

    Returns:
        Nusselt number
    """
    num = (f / 8) * (Re - 1000) * Pr
    den = 1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1)
    return num / den


def Nu_DittusBoelter(Re: float, Pr: float) -> float:
    """
    Calculate the Nusselt number using the Dittus-Boelter correlation.

    Args:
        Re: Reynolds number
        Pr: Prandtl number

    Returns:
        Nusselt number
    """
    return 0.023 * Re ** (4 / 5) * Pr ** (0.4)


def f_Pethukov(Re: float, Pr: float) -> float:
    """
    Calculate the friction factor using the Pethukov correlation.

    Args:
        Re: Reynolds number
        Pr: Prandtl number

    Returns:
        Friction factor
    """
    return (0.79 * np.log(Re) - 1.64) ** (-2)


def get_h_from_Nu(Nu: float, k: float, D: float) -> float:
    """
    Calculate the heat transfer coefficient from the Nusselt number.

    Args:
        Nu: Nusselt number
        k: Thermal conductivity
        D: Diameter

    Returns:
        Heat transfer coefficient
    """
    return Nu * k / D


def f_Haaland(Re: float, e_D: float) -> float:
    """
    Calculate the friction factor using the Haaland correlation.

    Args:
        Re: Reynolds number
        e_D: Relative roughness

    Returns:
        Friction factor
    """
    return (-1.8 * np.log10((e_D / 3.7) ** 1.11 + 6.9 / Re)) ** (-2)


def Schmidt(D: float, mu: float, rho: float) -> float:
    """
    Calculate the Schmidt number.

    Args:
        D: Diffusivity
        mu: Viscosity
        rho: Density

    Returns:
        Schmidt number
    """
    return mu / (rho * D)


def Sherwood(Sc: float, Re: float) -> float:
    """
    Calculate the Sherwood number.

    Args:
        Sc: Schmidt number
        Re: Reynolds number

    Returns:
        Sherwood number
    """
    return 0.0096 * Re**0.913 * Sc**0.346


def Sherwood_HT_analogy(Re: float, Sc: float) -> float:
    """
    Calculate the Sherwood number using the heat transfer analogy.

    Args:
        Re: Reynolds number
        Sc: Schmidt number

    Returns:
        Sherwood number
    """
    return 0.023 * Re**0.8 * Sc**0.4


def get_k_from_Sh(Sh: float, L: float, D: float) -> float:
    """
    Calculate the mass transfer coefficient from the Sherwood number.

    Args:
        Sh: Sherwood number
        L: Length
        D: Diameter

    Returns:
        Mass transfer coefficient
    """
    return Sh * D / L


def Re(rho: float, u: float, L: float, mu: float) -> float:
    """
    Calculate the Reynolds number.

    Args:
        rho: Density
        u: Velocity
        L: Length
        mu: Viscosity

    Returns:
        Reynolds number
    """
    return rho * u * L / mu


def Pr(c_p: float, mu: float, k: float) -> float:
    """
    Calculate the Prandtl number.

    Args:
        c_p: Specific heat capacity
        mu: Viscosity
        k: Thermal conductivity

    Returns:
        Prandtl number
    """
    return mu * c_p / k


def Sherwood_bubbles(Sc: float, Re: float) -> float:
    """
    Calculate the Sherwood number for bubbles.

    Args:
        Sc: Schmidt number
        Re: Reynolds number

    Returns:
        Sherwood number
    """
    return 0.089 * Re**0.69 * Sc**0.33


def get_length_HX(deltaTML: float, d_hyd: float, U: float, Q: float) -> float:
    """
    Calculate the length of the heat exchanger.

    Args:
        deltaTML: Log mean temperature difference
        d_hyd: Hydraulic diameter
        U: Overall heat transfer coefficient
        Q: Heat transfer rate

    Returns:
        Length of the heat exchanger
    """
    L = Q / (U * np.pi * d_hyd * deltaTML)
    return L


def get_deltaTML(T_in_hot: float, T_out_hot: float, T_in_cold: float, T_out_cold: float) -> float:
    """
    Calculate the log mean temperature difference.

    Args:
        T_in_hot: Hot inlet temperature
        T_out_hot: Hot outlet temperature
        T_in_cold: Cold inlet temperature
        T_out_cold: Cold outlet temperature

    Returns:
        Log mean temperature difference
    """
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    return delta_T_lm

