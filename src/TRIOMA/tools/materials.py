from TRIOMA.tools.component_tools import FluidMaterial, SolidMaterial
import numpy as np

N_A = 6.022e23
k_b = 8.617e-5
R_const = 8.314


# Define the fluid material
def Flibe(T):

    def density(T):
        return 2413 - 0.488 * T

    def viscosity(T):
        return 1.16e-4 * np.exp(3755 / T)

    def H_diff(T):
        R_const = 8.314
        return 9.3e-7 * np.exp(-42e3 / (R_const * T))

    def k():
        # thermal conductivity W/m/K
        return 1.1

    def cp():
        # soecific heat J/kg/K
        return 2386

    def k_H(T):
        # Henry's constant mol/m^3/Pa
        return 4.54e-4

    # def k_H_low(T):
    #     # Henry's constant mol/m^3/Pa Malinauskas
    #     R_const = 8.314
    #     return 2.1e-6 * np.exp(-29e3 / (R_const * T))

    # def k_H_high(T):
    #     # Henry's constant mol/m^3/Pa Calderoni
    #     R_const = 8.314
    #     return 7.9e-2 * np.exp(-35e3 / (R_const * T))

    # def k_H_mid(T):
    #     return 0.06 / T * np.exp(-(11 * (293 - 0.12 * T) / T))

    Flibe = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_H(T),
        MS=True,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )

    return Flibe


def Sodium(T):
    def density(T):
        return 219 + 275.32 * (1 - T / 2504.7) + 511.58 * (1 - T / 2503.7) ** 0.5

    def viscosity(T):
        return np.exp(-6.4406 - 0.3958 * np.log(T) + 556.835 / T)

    def H_diff(T):
        R_const = 8.314
        return 2e-5 * np.exp(-49053 / (R_const * T))

    def k():
        # thermal conductivity W/m/K
        return 1  # TODO

    def cp():
        # soecific heat J/kg/K
        return 1  # TODO

    def k_S(T):
        # Henry's constant mol/m^3/Pa
        return 10 ** (0.86 - 122 / T)

    Sodium = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_S(T),
        MS=False,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )
    return Sodium


def LiPb(T):
    def density(T):
        return 9659.8  # TODO

    def viscosity(T):
        return 1  # TODO

    def H_diff(T):
        return 1  # TODO

    def k():
        # thermal conductivity W/m/K
        return 1  # TODO

    def cp():
        # soecific heat J/kg/K
        return 1  # TODO

    def k_S(T):
        # Henry's constant mol/m^3/Pa
        MM_LiPb = 180  # TODO
        return 4.7e-7 * np.exp(-9e3 / (R_const * T) * density(T) / MM_LiPb)  # TODO

    LiPb = FluidMaterial(
        T,
        D=H_diff(T),
        Solubility=k_S(T),
        MS=False,
        mu=viscosity(T),
        rho=density(T),
        k=k(),
        cp=cp(),
    )
    return LiPb


def Steel(T):
    # Define the solid material

    def H_diff(T):
        D_met = 5.81e-7 * np.exp(-66.3e3 / (R_const * T))
        return D_met

    def K_S(T):
        return 1

    Steel = SolidMaterial(T=T, D=H_diff(T), K_S=K_S(T))

    return Steel
