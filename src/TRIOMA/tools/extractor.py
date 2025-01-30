import numpy
import TRIOMA.tools.correlations as cor
import scipy.integrate as integrate


def calculate_gas_velocity(G_gas, p_t, T, R):
    """__summary__
    Args:
        G_gas (float): Gas flowrate in normal conditions
        p_t (float): Total pressure
        T (float): Temperature
        R (float): Radius
    Returns:
        float: Gas velocity


    """
    Area = numpy.pi * R**2
    p_atm = 101325
    u_g = G_gas / Area / p_t * p_atm * T / 288.15
    return u_g


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


def length_extractor_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, kla, c_max=0):
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
    integral = NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, c_max)
    Z = u_l / kla * integral
    return Z


def NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, c_max=0):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity

    R_const = 8.314
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    R_g = 2 * u_g / u_l  ## gas on liquid ratio
    c_in = pl_in**0.5 * K_S  # inlet concentration in liquid
    c_out = pl_out**0.5 * K_S  # outlet concentration in liquid
    c_in_gas = pg_in / R_const / T

    def toint(c):
        value = 1 / (
            c - K_S * ((R_const * T / R_g) * (c - c_out + c_in_gas * R_g)) ** 0.5
        )
        return value

    c_g_max = pl_in / R_const / T  # maximum concentration in gas
    c_out_max = max(
        c_in - R_g * (c_g_max - pg_in / R_const / T),  # maximum gas stripping
        c_max,  # given input from equation
        (pg_in) ** 0.5 * K_S,  ## if liquid is in equilibrium with gas at outlet
    )

    integral = integrate.quad(toint, c_out, c_in, points=c_out_max, maxp1=1e3)
    if integral[0] < 0:
        ## for debugging reasons as for now
        print("Warning: negative value of the integral", integral)
        print(
            "c out is "
            + str(c_out)
            + " c in is "
            + str(c_in)
            + " c in gas is "
            + str(c_in_gas)
            + " p in liquid is "
            + str(pl_in)
            + " p in gas is "
            + str(pg_in)
            + "solubility is "
            + str(K_S)
            + "rg is "
            + str(R_g)
        )
        print(" c out max is " + str(c_out_max))
    return integral[0]


def NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, c_max=0):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    # but for molten salts by changing the evolution of c star following the Henry's law
    """
    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    c_in = pl_in * K_H
    c_out = pl_out * K_H
    R_const = 8.314
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    c_in_gas = pg_in / R_const / T
    R_g = u_g / u_l  ## gas on liquid ratio

    def toint(c):
        return 1 / (
            c - K_H * (u_l / u_g * R_const * T) * (c - c_out + c_in_gas * u_g / u_l)
        )

    c_g_max = pl_in / R_const / T  # maximum concentration in gas
    c_out_max = max(
        c_in - R_g * (c_g_max - pg_in / R_const / T),  # maximum gas stripping
        c_max,  # given input from equation
        pg_in * K_H,  ## if liquid is in equilibrium with gas at outlet
    )

    integral = integrate.quad(toint, c_out, c_in, points=c_out_max, maxp1=1e3)
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


def length_extractor_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, kla, c_max=0):
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
    integral = NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, c_max=c_max)

    Z = u_l / kla * integral
    return Z


from scipy.optimize import minimize, root


def get_c_out_GLC_lm(Z, R, G_l, G_gas, pl_in, T, p_t, K_S, pg_in, kla):
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
    u_l = G_l / (numpy.pi * R**2)
    c_in = pl_in**0.5 * K_S
    c_out_max_reaction = c_in - kla * Z / u_l * (
        c_in - pg_in**0.5 * K_S
    )  ## concentration at the outlet assuming maximum reaction rate possible
    R_const = 8.314
    Area = numpy.pi * R**2
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    c_g_max = pl_in / R_const / T  ## maximum gas concentration according with liquid
    c_out_max = max(
        c_in
        - u_g
        / u_l
        * 2
        * (
            c_g_max - pg_in / R_const / T
        ),  ## liquid concentration if gas strips as much as possible and gets into eq with liquid
        (pg_in) ** 0.5
        * K_S,  ## liquid concentration if it gets in equilibrium with gas
        c_out_max_reaction,  ## liquid concentration if reaction rate is at the maximum
    )

    def lenght_residual(c_out):
        pl_out_2 = c_out**2 / K_S**2
        z_guess = length_extractor_lm(
            R,
            G_l,
            G_gas,
            pl_in,
            pl_out_2,
            T,
            p_t,
            K_S,
            pg_in,
            kla,
            c_max=c_out_max,
        )
        return abs(float(Z - z_guess) ** 2)

    c_out = minimize(
        lenght_residual,
        c_out_max + (c_in - c_out_max) / 2,
        method="Powell",
        bounds=[(float(c_out_max), float(c_in))],
        tol=1e-20,
        options={"maxiter": 1e8, "xatol": 1e-20, "fatol": 1e-20},
    ).x[0]
    eff = 1 - c_out / c_in
    L_cout = length_extractor_lm(
        R, G_l, G_gas, pl_in, c_out**2 / K_S**2, T, p_t, K_S, pg_in, kla
    )
    if abs(L_cout - Z) > 1e-3:
        print(
            "Warning!: guessed length is not equal to the height. Double check your result",
            L_cout,
            Z,
        )
    if abs(c_out - c_out_max) < 1e-8:
        print("The sweep gas saturated")
        if L_cout < Z:
            print(" Longer column would not increment the extraction efficiency")
            eff = 1 - c_out / c_in
    return c_out, eff


def get_c_out_GLC_ms(Z, R, G_l, G_gas, pl_in, T, p_t, K_H, pg_in, kla):
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
    u_l = G_l / (numpy.pi * R**2)
    c_in = pl_in * K_H
    c_out_max_reaction = c_in - kla * Z / u_l * (
        c_in - pg_in * K_H
    )  ## concentration at the outlet assuming maximum reaction rate possible
    R_const = 8.314
    Area = numpy.pi * R**2
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    c_g_max = pl_in / R_const / T  ## maximum gas concentration according with liquid
    c_out_max = max(
        c_in
        - u_g
        / u_l
        * (
            c_g_max - pg_in / R_const / T
        ),  ## liquid concentration if gas strips as much as possible and gets into eq with liquid
        (pg_in) * K_H,  ## liquid concentration if it gets in equilibrium with gas
        c_out_max_reaction,  ## liquid concentration if reaction rate is at the maximum
    )

    def lenght_residual(c_out):
        pl_out_2 = c_out / K_H
        z_guess = length_extractor_ms(
            R,
            G_l,
            G_gas,
            pl_in,
            pl_out_2,
            T,
            p_t,
            K_H,
            pg_in,
            kla,
            c_max=c_out_max,
        )
        return abs(float(Z - z_guess) ** 2)

    c_out = minimize(
        lenght_residual,
        c_out_max + (c_in - c_out_max) / 2,
        method="Powell",
        bounds=[(float(c_out_max), float(c_in))],
        tol=1e-20,
        options={"maxiter": 1e8, "xatol": 1e-20, "fatol": 1e-20},
    ).x[0]
    eff = 1 - c_out / c_in
    L_cout = length_extractor_lm(
        R, G_l, G_gas, pl_in, c_out / K_H, T, p_t, K_H, pg_in, kla
    )
    if abs(L_cout - Z) > 1e-3:
        print(
            "Warning!: guessed length is not equal to the height. Double check your result",
            L_cout,
            Z,
        )
    if abs(c_out - c_out_max) < 1e-8:
        print("The sweep gas saturated")
        if L_cout < Z:
            print(" Longer column would not increment the extraction efficiency")
            eff = 1 - c_out / c_in
    return c_out, eff


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
