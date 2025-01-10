import numpy
import tools.correlations as cor
import scipy.integrate as integrate


def extractor_lm(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in):
    """_summary_
    Args:
        Z (float): Height
        R (float): Radius
        G_l (float): Liquid flowrate
        G_gas (float): Gas flowrate
        pl_in (float): T pressure inlet
        pl_out (float): T pressure outlet
        T (float): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2  ## Area of the column
    B_l = (G_l) / Area * 3600  ## Liquid holdup
    u_l = G_l / Area  # Liquid velocity
    integral = NTU_lm(
        R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in
    )  ## number of transfer units,
    kla_c = u_l / Z * integral
    return [B_l, kla_c]


def length_extractor_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, kla):
    """_summary_
    Args:
        Z (float): Height
        R (float): Radius
        G_l (float): Liquid flowrate
        G_gas (float): Gas flowrate
        pl_in (float): T pressure inlet
        pl_out (float): T pressure outlet
        T (float): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    integral = NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in)
    Z = u_l / kla * integral
    if Z < 0:
        print("Warning: the system is not feasable or the integral is not accurate")
        return None
    return Z


def NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    R_g = 2 * G_gas / G_l
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15
    c_in = pl_in**0.5 * K_S
    c_out = pl_out**0.5 * K_S
    c_in_gas = pg_in / R_const / T

    def toint(c):
        return 1 / (
            c
            - K_S
            * (u_l / 2 / u_g * R_const * T) ** 0.5
            * (c - c_out + c_in_gas * 2 * u_g / u_l) ** 0.5
        )

    integral = integrate.quad(toint, c_out, c_in)
    return integral[0]


def NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    # but for molten salts by changing the evolution of c star following the Henry's law
    """
    Area = numpy.pi * R**2
    B_l = (G_l) / Area * 3600
    u_l = G_l / Area  # Liquid velocity
    c_in = pl_in * K_H
    c_out = pl_out * K_H
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15
    c_in_gas = pg_in / R_const / T

    def toint(c):
        return 1 / (
            c - K_H * (u_l / u_g * R_const * T) * (c - c_out + c_in_gas * u_g / u_l)
        )

    integral = integrate.fixed_quad(toint, c_out, c_in)
    return integral[0]


def extractor_ms(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in):
    """_summary_
    Args:
        Z (float): Height
        R (float): Radius
        G_l (float): Liquid flowrate
        G_gas (float): Gas flowrate
        pl_in (float): T pressure inlet
        pl_out (float): T pressure outlet
        T (float): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    B_l = (G_l) / Area * 3600
    u_l = G_l / Area  # Liquid velocity
    integral = NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in)
    kla_c = u_l / Z * integral
    return [B_l, kla_c]


def length_extractor_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, kla):
    """_summary_
    Args:
        R (float): Radius
        G_l (float): Liquid flowrate
        G_gas (float): Gas flowrate
        pl_in (float): T pressure inlet
        pl_out (float): T pressure outlet
        T (float): Temperature
        p_t(float): pressure Pa of the column
        K_S(float): Sievert's constant
    Returns:
        B_l (float):liquid load
        k_la(float): mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    integral = NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in)

    Z = u_l / kla * integral
    if Z < 0:
        print("Warning: the system is not feasable or the integral is not accurate")
        return None
    return Z


def pack_corr(a, d, D, eta, v):
    k_l = (
        0.0051
        * (v / eta / a) ** (2 / 3)
        * (D / eta) ** 0.5
        * (a * d) ** 0.4
        * (1 / eta / 9.81) ** (-1 / 3)
    )
    return k_l


def corr_packed(Re, Sc, d, rho_L, mu_L, L, D):
    """_summary_
    Args
        Re (float): Reynolds
        Sc (float): Schmidt
        d (float): ring diameter
        rho_L (float): Liquid density
        mu_L (float): viscosity
        D diffusion coeff
        L= characteristic length
    Returns:
        float: _description_
        Warning: verification of this must be done
    """
    beta = 0.32  # Raschig rings 0.25
    g = 9.81
    Sh = beta * Re**0.59 * Sc**0.5 * (d**3 * g * rho_L**2 / mu_L**2) ** 0.17
    return cor.get_k_from_Sh(Sh, L, D)
