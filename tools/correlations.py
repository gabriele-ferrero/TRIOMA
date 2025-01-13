import numpy as np


def Nu_SiederTate(Re, Pr, mu_w, mu_c):
    """
    Calculate the Nusselt number using the Sieder-Tate correlation.

    Parameters:
    - Re (float): Reynolds number
    - Pr (float): Prandtl number
    - mu_w (float): viscosity of the wall fluid
    - mu_c (float): viscosity of the core fluid

    Returns:
    - Nu (float): Nusselt number
    """
    return 0.027 * Re ** (4 / 5) * Pr ** (0.3) * (mu_c / mu_w) ** (0.14)


def Nu_Gnielinsky(Re, Pr, f):
    """Calculate the Nusselt number using the Gnielinsky correlation.
    Parameters:
    - Re (float): Reynolds number
    - Pr (float): Prandtl number
    - f (float): friction factor
    Returns:
    - Nu (float): Nusselt number
    """
    num = (f / 8) * (Re - 1000) * Pr
    den = 1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1)
    return num / den


def Nu_DittusBoelter(Re, Pr):
    """Calculate the Nusselt number using the Dittus-Boelter correlation.
    Parameters:
    - Re (float): Reynolds number
    - Pr (float): Prandtl number
    Returns:
    - Nu (float): Nusselt number
    """
    return 0.023 * Re ** (4 / 5) * Pr ** (0.4)


def f_Pethukov(Re, Pr):
    """Calculate the friction factor using the Pethukov correlation.
    Parameters:
    - Re (float): Reynolds number
    - Pr (float): Prandtl number
    Returns:
    - f (float): friction factor
    """
    return (0.79 * np.log(Re) - 1.64) ** (-2)


def get_h_from_Nu(Nu, k, D):
    """Calculate the heat transfer coefficient from the Nusselt number.
    Parameters:
    - Nu (float): Nusselt number
    - k (float): thermal conductivity
    - D (float): diameter
    Returns:
    - h (float): heat transfer coefficient
    """
    return Nu * k / D


def f_Haaland(Re, e_D):
    """Calculate the friction factor using the Haaland correlation.
    Parameters:
    - Re (float): Reynolds number
    - e_D (float): relative roughness
    Returns:
    - f (float): friction factor
    """
    return (-1.8 * np.log10((e_D / 3.7) ** 1.11 + 6.9 / Re)) ** (-2)


def Schmidt(D, mu, rho):
    """Calculate the Schmidt number.
    Parameters:
    - D (float): diffusivity
    - mu (float): viscosity
    - rho (float): density
    Returns:
    - Sc (float): Schmidt number
    """
    return mu / (rho * D)


def Sherwood(Sc, Re):
    """Calculate the Sherwood number.
    Parameters:
    - Sc (float): Schmidt number
    - Re (float): Reynolds number
    Returns:
    - Sh (float): Sherwood number
    """

    return 0.0096 * Re**0.913 * Sc**0.346


def Sherwood_HT_analogy(Re, Sc):
    """Calculate the Sherwood number using the heat transfer analogy.
    Parameters:
    - Re (float): Reynolds number
    - Sc (float): Schmidt number
    Returns:
    - Sh (float): Sherwood number
    """
    return 0.023 * Re**0.8 * Sc**0.4


def get_k_from_Sh(Sh, L, D):
    """Calculate the mass transfer coefficient from the Sherwood number.
    Parameters:
    - Sh (float): Sherwood number
    - L (float): length
    - D (float): diameter
    Returns:
    - k (float): mass transfer coefficient
    """
    return Sh * D / L


def Re(rho, u, L, mu):
    """Calculate the Reynolds number.
    Parameters:
    - rho (float): density
    - u (float): velocity
    - L (float): length
    - mu (float): viscosity
    Returns:
    - Re (float): Reynolds number
    """
    return rho * u * L / mu


def Pr(c_p, mu, k):
    """Calculate the Prandtl number.
    Parameters:
    - c_p (float): specific heat capacity
    - mu (float): viscosity
    - k (float): thermal conductivity
    Returns:
    - Pr (float): Prandtl number
    """
    return mu * c_p / k


def Sherwood_bubbles(Sc, Re):
    """Calculate the Sherwood number for bubbles.
    Parameters:
    - Sc (float): Schmidt number
    - Re (float): Reynolds number
    Returns:
    - Sh (float): Sherwood number
    """
    return 0.089 * Re**0.69 * Sc**0.33  ## Humrickhouse MS report


def get_length_HX(deltaTML, d_hyd, U, Q):
    """Calculate the length of the heat exchanger.
    Parameters:
    - deltaTML (float): log mean temperature difference
    - d_hyd (float): hydraulic diameter
    - U (float): overall heat transfer coefficient
    - Q (float): heat transfer rate
    Returns:
    - L (float): length of the heat exchanger
    """
    L = Q / (U * np.pi * d_hyd * deltaTML)
    return L


def get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold):
    """Calculate the log mean temperature difference.
    Parameters:
    - T_in_hot (float): hot inlet temperature
    - T_out_hot (float): hot outlet temperature
    - T_in_cold (float): cold inlet temperature
    - T_out_cold (float): cold outlet temperature
    Returns:
    - deltaTML (float): log mean temperature difference
    """
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    return delta_T_lm
