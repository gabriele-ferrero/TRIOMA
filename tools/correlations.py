import numpy as np
import matplotlib.pyplot as plt


def Nu_SiederTate(Re, Pr, mu_w, mu_c):
    return 0.027 * Re ** (4 / 5) * Pr ** (0.3) * (mu_c / mu_w) ** (0.14)


def Nu_Gnielinsky(Re, Pr, f):
    num = (f / 8) * (Re - 1000) * Pr
    den = 1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1)
    return num / den


def Nu_DittusBoelter(Re, Pr):
    return 0.023 * Re ** (4 / 5) * Pr ** (0.4)


def f_Pethukov(Re, Pr):
    return (0.79 * np.log(Re) - 1.64) ** (-2)


def get_h_from_Nu(Nu, k, D):
    return Nu * k / D


def f_Haaland(Re, e_D):
    return (-1.8 * np.log10((e_D / 3.7) ** 1.11 + 6.9 / Re)) ** (-2)


def Schmidt(D, mu, rho):
    return mu / (rho * D)


def Sherwood(Sc, Re):
    return 0.0096 * Re**0.913 * Sc**0.346


def Sherwood_HT_analogy(Re, Sc):
    return 0.023 * Re**0.8 * Sc**0.4


def get_k_from_Sh(Sh, L, D):
    return Sh * D / L


def Re(rho, u, L, mu):
    return rho * u * L / mu


def Pr(c_p, mu, k):
    return mu * c_p / k


def Sherwood_bubbles(Sc, Re):
    return 0.089 * Re**0.69 * Sc**0.33  ## Humrickhouse MS report


def get_length_HX(deltaTML, d_hyd, U, Q):
    L = Q / (U * np.pi * d_hyd * deltaTML)
    return L


def get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold):
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    return delta_T_lm
