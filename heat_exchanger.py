import numpy as np
from tools.component_tools import Component
from tools.component_tools import Fluid
from tools.component_tools import Membrane
from tools.component_tools import GLC_Gas
from tools.component_tools import GLC
from tools.materials import Flibe
import tools.materials as materials
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy.optimize import fsolve

# Define known parameters
U = 500  # Overall heat transfer coefficient in W/m^2.K
A = 50  # Area in m^2
T_in_hot = 150  # Hot fluid inlet temperature in Celsius
T_out_cold = 30  # Cold fluid outlet temperature in Celsius
m_dot_cold = 2  # Mass flow rate of cold fluid in kg/s
c_p_cold = 4184  # Specific heat capacity of cold fluid in J/kg.K


# def equations(p):
#     T_out_hot, T_in_cold = p
#     Q = (
#         m_dot_cold * c_p_cold * (T_out_cold - T_in_cold)
#     )  # Heat transfer from cold fluid
#     delta_T1 = T_in_hot - T_out_cold
#     delta_T2 = T_out_hot - T_in_cold
#     delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
#     Q_calculated = U * A * delta_T_lm
#     return (Q - Q_calculated, T_in_hot - T_out_hot - (T_in_cold - T_out_cold))


# # Initial guesses for T_out_hot and T_in_cold
# initial_guesses = [120, 20]
# solution = fsolve(equations, initial_guesses)
# print(f"T_out_hot: {solution[0]} °C, T_in_cold: {solution[1]} °C")
c0 = 4.1e-4
T = 800
d_hyd = 1e-3

mat = materials.Flibe(T)

flibe = Fluid(
    T=T,
    Solubility=mat.Solubility,
    MS=True,
    D=mat.D,
    d_Hyd=d_hyd,
    mu=mat.mu,
    rho=mat.rho,
    U0=1,
    c0=c0,
)
Steel = Membrane(T=T, D=1e-9, thick=1e-3, K_S=1e-3, k_d=1e6, k_r=1e6)
PAV = Component(c_in=c0, eff=0.5, fluid=flibe, membrane=Steel)
U = PAV.get_global_HX_coeff()


def get_length(deltaTML, d_hyd, U, Q):
    L = Q / (U * np.pi * d_hyd * deltaTML)
    return L


def get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold):
    delta_T1 = T_in_hot - T_out_cold
    delta_T2 = T_out_hot - T_in_cold
    delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    return delta_T_lm


L = get_length(get_deltaTML(150, 120, 20, 30), d_hyd, PAV.U, 1e6)
