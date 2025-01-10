import numpy
import tools.correlations as cor
import scipy.integrate as integrate


def extractor_lm(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in):
    """_summary_
    Args:
        Z (_type_): Height
        R (_type_): Radius
        G_l (_type_): Liquid flowrate
        G_gas (_type_): Gas flowrate
        pl_in (_type_): T pressure inlet
        pl_out (_type_): T pressure outlet
        T (_type_): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    B_l = (G_l) / Area * 3600
    u_l = G_l / Area  # Liquid velocity
    c_in = pl_in**0.5 * K_S
    c_out = pl_out**0.5 * K_S
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15

    def toint(c):
        return 1 / (c - K_S * (u_l / 2 / u_g * R_const * T) ** 0.5 * (c - c_out) ** 0.5)

    integral = integrate.quad(toint, c_out, c_in)
    kla_c = u_l / Z * integral[0]
    return [B_l, kla_c]


def length_extractor_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, kla):
    """_summary_
    Args:
        Z (_type_): Height
        R (_type_): Radius
        G_l (_type_): Liquid flowrate
        G_gas (_type_): Gas flowrate
        pl_in (_type_): T pressure inlet
        pl_out (_type_): T pressure outlet
        T (_type_): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    R_g = 2 * G_gas / G_l
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15
    c_in = pl_in**0.5 * K_S
    c_out = pl_out**0.5 * K_S

    def toint(c):
        return 1 / (c - K_S * (u_l / 2 / u_g * R_const * T) ** 0.5 * (c - c_out) ** 0.5)

    integral = integrate.quad(toint, c_out, c_in)
    Z = u_l / kla * integral[0]
    return Z


def extractor_ms(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in):
    """_summary_
    Args:
        Z (_type_): Height
        R (_type_): Radius
        G_l (_type_): Liquid flowrate
        G_gas (_type_): Gas flowrate
        pl_in (_type_): T pressure inlet
        pl_out (_type_): T pressure outlet
        T (_type_): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    B_l = (G_l) / Area * 3600
    u_l = G_l / Area  # Liquid velocity
    c_in = pl_in * K_H
    c_out = pl_out * K_H
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15

    def toint(c):
        return 1 / (c - K_H * (u_l / u_g * R_const * T) * (c - c_out))

    integral = integrate.fixed_quad(toint, c_out, c_in)
    print("integral is", integral)
    kla_c = u_l / Z * integral[0]
    return [B_l, kla_c]


def length_extractor_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, kla):
    """_summary_
    Args:
        Z (_type_): Height
        R (_type_): Radius
        G_l (_type_): Liquid flowrate
        G_gas (_type_): Gas flowrate
        pl_in (_type_): T pressure inlet
        pl_out (_type_): T pressure outlet
        T (_type_): Temperature
        p_t pressure Pa of the column
        K_S Sievert's constant
    Returns:
        B_l liquid load
        k_la mass transfer coefficient in packed column
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    R_g = 2 * G_gas / G_l
    R_const = 8.314
    u_g = G_gas / Area / p_t * 1e5 * T / 288.15
    c_in = pl_in * K_H
    c_out = pl_out * K_H

    def toint(c):
        return 1 / (c - K_H * (u_l / u_g * R_const * T) * (c - c_out))

    integral = integrate.fixed_quad(toint, c_out, c_in)
    Z = u_l / kla * integral[0]
    if Z < 0:
        print("Warning: the system is not feasable")
        return 0
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
        Re (_type_): Reynolds
        Sc (_type_): Schmidt
        d (_type_): ring diameter
        rho_L (_type_): Liquid density
        mu_L (_type_): viscosity
        D diffusion coeff
        L= characteristic length
    Returns:
        _type_: _description_
        Warning: verification of this must be done
    """
    beta = 0.32  # Raschig rings 0.25
    g = 9.81
    Sh = beta * Re**0.59 * Sc**0.5 * (d**3 * g * rho_L**2 / mu_L**2) ** 0.17
    return cor.get_k_from_Sh(Sh, L, D)


# def length_extractor_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, kla):
#     """_summary_
#     Args:
#         Z (_type_): Height
#         R (_type_): Radius
#         G_l (_type_): Liquid flowrate
#         G_gas (_type_): Gas flowrate
#         pl_in (_type_): T pressure inlet
#         pl_out (_type_): T pressure outlet
#         T (_type_): Temperature
#         p_t pressure Pa of the column
#         K_S Sievert's constant
#     Returns:
#         B_l liquid load
#         k_la mass transfer coefficient in packed column
#     """
#     Area = numpy.pi * R**2
#     B_l = (G_l) / Area * 3600
#     u_l = G_l / Area  # Liquid velocity
#     R_g = G_gas / G_l
#     A_p = (K_S * p_t * 0.0224 / R_g) ** 0.5 * (
#         pl_in**0.5 - pl_out**0.5 + pg_in * R_g / (K_S * p_t * 0.0224)
#     ) ** 0.5

#     Z = u_l / kla * numpy.log((pl_in**0.5 - A_p) / (pl_out**0.5 - A_p))
#     return Z
