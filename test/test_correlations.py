import sys
import os

sys.path.append(os.path.abspath("."))
from src.TRIOMA.tools.correlations import (
    get_deltaTML,
    Nu_SiederTate,
    Nu_Gnielinsky,
    Nu_DittusBoelter,
    f_Pethukov,
    f_Haaland,
    get_h_from_Nu,
    Schmidt,
    Re,
    Pr,
    get_k_from_Sh,
    Sherwood,
    Sherwood_HT_analogy,
    get_length_HX,
    Sherwood_bubbles,
)
import numpy as np


def test_get_deltaTML():
    T_in_hot = 100.0
    T_out_hot = 80.0
    T_in_cold = 40.0
    T_out_cold = 20.0

    expected_deltaTML = ((T_in_hot - T_out_cold) - (T_out_hot - T_in_cold)) / np.log(
        (T_in_hot - T_out_cold) / (T_out_hot - T_in_cold)
    )

    deltaTML = get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold)

    assert deltaTML == expected_deltaTML, "Test failed: Incorrect deltaTML value"


def test_Nu_SiederTate():
    Re = 10000.0
    Pr = 0.7
    mu_w = 0.001
    mu_c = 0.002

    expected_Nu = 0.027 * Re ** (4 / 5) * Pr ** (0.3) * (mu_c / mu_w) ** (0.14)

    Nu = Nu_SiederTate(Re, Pr, mu_w, mu_c)

    assert Nu == expected_Nu, "Test failed: Incorrect Nusselt number value"


def test_Nu_Gnielinsky():
    Re = 10000.0
    Pr = 0.7
    f = 0.02

    expected_Nu = (
        (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    )

    Nu = Nu_Gnielinsky(Re, Pr, f)

    assert Nu == expected_Nu, "Test failed: Incorrect Nusselt number value"


def test_Nu_DittusBoelter():
    Re = 10000.0
    Pr = 0.7

    expected_Nu = 0.023 * Re ** (4 / 5) * Pr ** (0.4)

    Nu = Nu_DittusBoelter(Re, Pr)

    assert Nu == expected_Nu, "Test failed: Incorrect Nusselt number value"


def test_f_Pethukov():
    Re = 10000.0
    Pr = 0.7

    expected_f = (0.79 * np.log(Re) - 1.64) ** (-2)

    f = f_Pethukov(Re, Pr)

    assert f == expected_f, "Test failed: Incorrect friction factor value"


def test_get_h_from_Nu():
    Nu = 100.0
    k = 0.5
    D = 0.1

    expected_h = Nu * k / D

    h = get_h_from_Nu(Nu, k, D)

    assert h == expected_h, "Test failed: Incorrect heat transfer coefficient value"


def test_f_Haaland():
    Re = 10000.0
    e_D = 0.001

    expected_f = (-1.8 * np.log10((e_D / 3.7) ** 1.11 + 6.9 / Re)) ** (-2)

    f = f_Haaland(Re, e_D)

    assert f == expected_f, "Test failed: Incorrect friction factor value"


def test_Schmidt():
    D = 0.1
    mu = 0.001
    rho = 1000.0

    expected_Sc = mu / (rho * D)

    Sc = Schmidt(D, mu, rho)

    assert Sc == expected_Sc, "Test failed: Incorrect Schmidt number value"


def test_Sherwood():
    Sc = 0.7
    Re = 10000.0

    expected_Sh = 0.0096 * Re**0.913 * Sc**0.346

    Sh = Sherwood(Sc, Re)

    assert Sh == expected_Sh, "Test failed: Incorrect Sherwood number value"


def test_Sherwood_HT_analogy():
    Re = 10000.0
    Sc = 0.7

    expected_Sh = 0.023 * Re**0.8 * Sc**0.4

    Sh = Sherwood_HT_analogy(Re, Sc)

    assert Sh == expected_Sh, "Test failed: Incorrect Sherwood number value"


def test_get_k_from_Sh():
    Sh = 100.0
    L = 1.0
    D = 0.1

    expected_k = Sh * D / L

    k = get_k_from_Sh(Sh, L, D)

    assert k == expected_k, "Test failed: Incorrect mass transfer coefficient value"


def test_Re():
    rho = 1000.0
    u = 2.0
    L = 0.5
    mu = 0.001

    expected_Re = rho * u * L / mu

    Re_val = Re(rho, u, L, mu)

    assert Re_val == expected_Re, "Test failed: Incorrect Reynolds number value"


def test_Pr():
    c_p = 1000.0
    mu = 0.001
    k = 0.5

    expected_Pr = mu * c_p / k

    Pr_val = Pr(c_p, mu, k)

    assert Pr_val == expected_Pr, "Test failed: Incorrect Prandtl number value"


def test_get_length_HX():
    deltaTML = 20.0
    d_hyd = 0.1
    U = 100.0
    Q = 5000.0

    expected_L = Q / (U * np.pi * d_hyd * deltaTML)

    L = get_length_HX(deltaTML, d_hyd, U, Q)

    assert L == expected_L, "Test failed: Incorrect length of the heat exchanger value"


def test_Sherwood_bubbles():
    Sc = 0.7
    Re = 10000.0

    expected_Sh = 0.089 * Re**0.69 * Sc**0.33

    Sh = Sherwood_bubbles(Sc, Re)

    assert Sh == expected_Sh, "Test failed: Incorrect Sherwood number value"
