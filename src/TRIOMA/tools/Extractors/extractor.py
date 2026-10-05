import numpy
import TRIOMA.tools.correlations as cor
import scipy.integrate as integrate
from scipy.optimize import brentq


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


def _integrate_linear_denominator(c_out, c_in, slope, intercept):
    """
    Integrate 1 / (slope*c + intercept) from c_out to c_in.

    Returns infinity when the denominator vanishes at an endpoint,
    which corresponds to an equilibrium/saturation limit.
    """

    d_out = slope * c_out + intercept
    d_in = slope * c_in + intercept

    tolerance = 1.0e-14 * max(
        1.0,
        abs(d_out),
        abs(d_in),
        abs(c_out),
        abs(c_in),
    )

    # Endpoint singularity: the physical equilibrium limit.
    if abs(d_out) <= tolerance or abs(d_in) <= tolerance:
        return numpy.inf

    # A linear denominator changing sign has a pole inside the interval.
    if d_out * d_in < 0.0:
        raise ValueError(
            "The MS GLC NTU integrand has a singularity inside " "the integration interval."
        )

    # Constant denominator.
    if abs(slope) <= numpy.finfo(float).eps:
        if abs(intercept) <= tolerance:
            return numpy.inf
        return (c_in - c_out) / intercept

    integral_value = numpy.log(abs(d_in / d_out)) / slope

    if numpy.isnan(integral_value):
        raise ValueError("The MS GLC NTU integral returned NaN.")

    return float(integral_value)


def NTU_lm(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S, pg_in, c_max=0):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    """
    # Convert array inputs to scalars
    pl_out = numpy.asarray(pl_out).item() if numpy.asarray(pl_out).ndim > 0 else float(pl_out)
    c_max = numpy.asarray(c_max).item() if numpy.asarray(c_max).ndim > 0 else float(c_max)

    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity

    R_const = 8.314
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    R_g = 2 * u_g / u_l  ## gas on liquid ratio
    c_in = float(pl_in**0.5 * K_S)  # inlet concentration in liquid
    c_out = float(pl_out**0.5 * K_S)  # outlet concentration in liquid
    c_in_gas = float(pg_in / R_const / T)

    def toint(c):
        value = 1 / (c - K_S * ((R_const * T / R_g) * (c - c_out + c_in_gas * R_g)) ** 0.5)
        return value

    c_g_max = float(pl_in / R_const / T)  # maximum concentration in gas
    c_out_max = float(
        max(
            c_in - R_g * (c_g_max - pg_in / R_const / T),  # maximum gas stripping
            c_max,  # given input from equation
            (pg_in) ** 0.5 * K_S,  ## if liquid is in equilibrium with gas at outlet
        )
    )

    integral = integrate.quad(toint, c_out, c_in, points=c_out_max, maxp1=1e3)
    if integral[0] < 0 or not numpy.isfinite(integral[0]):
        # quad cannot resolve the endpoint singularity when c_out approaches its
        # equilibrium limit (e.g. pg_in=0 -> c_out_max=0): the true NTU is divergent
        # (infinite length needed to reach that limit exactly), not negative.
        return numpy.inf
    return integral[0]


def NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, c_max=0):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    # but for molten salts by changing the evolution of c star following the Henry's law
    """
    # Convert array inputs to scalars
    pl_out = numpy.asarray(pl_out).item() if numpy.asarray(pl_out).ndim > 0 else float(pl_out)
    c_max = numpy.asarray(c_max).item() if numpy.asarray(c_max).ndim > 0 else float(c_max)

    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    c_in = float(pl_in * K_H)
    c_out = float(pl_out * K_H)
    R_const = 8.314
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    c_in_gas = float(pg_in / R_const / T)
    R_g = u_g / u_l  ## gas on liquid ratio

    def toint(c):
        return 1 / (c - K_H * (u_l / u_g * R_const * T) * (c - c_out + c_in_gas * u_g / u_l))

    c_g_max = float(pl_in / R_const / T)  # maximum concentration in gas
    c_out_max = float(
        max(
            c_in - R_g * (c_g_max - pg_in / R_const / T),  # maximum gas stripping
            c_max,  # given input from equation
            pg_in * K_H,  ## if liquid is in equilibrium with gas at outlet
        )
    )

    # integral = integrate.quad(toint, c_out, c_in, points=c_out_max, maxp1=1e3)
    integral = integrate.fixed_quad(toint, c_out, c_in)
    return integral[0]


def NTU_ms(R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H, pg_in, c_max=0):
    """
    solving integral equation from (5) of "The engineering sizing of the packed desorption column of hydrogen
    # isotopes from Pb–17Li eutectic alloy. A rate based model using
    # experimental mass transfer coefficients from a Melodie loop""
    # but for molten salts by changing the evolution of c star following the Henry's law
    """
    # Convert array inputs to scalars
    pl_out = numpy.asarray(pl_out).item() if numpy.asarray(pl_out).ndim > 0 else float(pl_out)
    c_max = numpy.asarray(c_max).item() if numpy.asarray(c_max).ndim > 0 else float(c_max)

    Area = numpy.pi * R**2
    u_l = G_l / Area  # Liquid velocity
    c_in = float(pl_in * K_H)
    c_out = float(pl_out * K_H)
    R_const = 8.314
    u_g = calculate_gas_velocity(G_gas=G_gas, p_t=p_t, T=T, R=R)
    c_in_gas = float(pg_in / R_const / T)
    R_g = u_g / u_l  ## gas on liquid ratio

    def toint(c):
        return 1 / (c - K_H * (u_l / u_g * R_const * T) * (c - c_out + c_in_gas * u_g / u_l))

    c_g_max = float(pl_in / R_const / T)  # maximum concentration in gas
    c_out_max = float(
        max(
            c_in - R_g * (c_g_max - pg_in / R_const / T),  # maximum gas stripping
            c_max,  # given input from equation
            pg_in * K_H,  ## if liquid is in equilibrium with gas at outlet
        )
    )

    # Physical outlet-concentration checks.
    concentration_tolerance = 1.0e-12 * max(1.0, abs(c_in))

    if c_out_max > c_in + concentration_tolerance:
        raise ValueError(
            "No physical MS GLC extraction solution exists: "
            f"c_out_max={c_out_max:.6e} is greater than "
            f"c_in={c_in:.6e}."
        )

    if c_out > c_in + concentration_tolerance:
        raise ValueError(
            f"c_out={c_out:.6e} cannot be greater than "
            f"c_in={c_in:.6e} for an extraction calculation."
        )

    if c_out < c_out_max - concentration_tolerance:
        raise ValueError(
            "The requested MS outlet concentration is below the "
            "physical gas-saturation/equilibrium limit: "
            f"c_out={c_out:.6e}, c_out_min={c_out_max:.6e}."
        )

    # No concentration change means zero NTU.
    if c_out >= c_in - concentration_tolerance:
        return 0.0

    # The MS integrand has the form:
    #
    #     1 / (slope * c + intercept)
    #
    # where
    #
    #     c - K_H * (u_l / u_g * R_const * T)
    #         * (c - c_out + c_in_gas * u_g / u_l)
    #
    # is linear in c.

    A = K_H * (u_l / u_g) * R_const * T

    slope = 1.0 - A
    intercept = A * (c_out - c_in_gas * u_g / u_l)

    integral_value = _integrate_linear_denominator(
        c_out=c_out,
        c_in=c_in,
        slope=slope,
        intercept=intercept,
    )

    # Infinite NTU is expected at the physical saturation limit.
    # It is handled by get_c_out_GLC_ms().
    if numpy.isnan(integral_value):
        raise ValueError("The MS GLC NTU integral returned NaN.")
    if integral_value < 0.0:
        print(
            "Warning: negative MS GLC NTU integral. "
            f"NTU={integral_value:.6e}, "
            f"c_out={c_out:.6e}, "
            f"c_in={c_in:.6e}, "
            f"c_out_min={c_out_max:.6e}."
        )

    return integral_value


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


def get_c_out_GLC_lm(
    Z,
    R,
    G_l,
    G_gas,
    pl_in,
    T,
    p_t,
    K_S,
    pg_in,
    kla,
):
    """
    Calculate the LM GLC outlet concentration and extraction efficiency.

    Parameters
    ----------
    Z : float
        GLC height [m].
    R : float
        GLC radius [m].
    G_l : float
        Liquid volumetric flow rate [m^3/s].
    G_gas : float
        Gas flow rate at reference conditions [m^3/s].
    pl_in : float
        Liquid inlet isotope partial pressure [Pa].
    T : float
        GLC temperature [K].
    p_t : float
        Total gas pressure [Pa].
    K_S : float
        Sievert solubility constant.
    pg_in : float
        Gas inlet isotope partial pressure [Pa].
    kla : float
        Volumetric liquid-side mass-transfer coefficient [1/s].

    Returns
    -------
    c_out : float
        Liquid outlet isotope concentration [mol/m^3].
    eff : float
        Liquid extraction efficiency.
    """

    if Z < 0.0:
        raise ValueError("The GLC height Z must be non-negative.")
    if R <= 0.0:
        raise ValueError("The GLC radius R must be positive.")
    if G_l <= 0.0:
        raise ValueError("The liquid flow rate G_l must be positive.")
    if G_gas <= 0.0:
        raise ValueError("The gas flow rate G_gas must be positive.")
    if kla < 0.0:
        raise ValueError("The mass-transfer coefficient kla cannot be negative.")

    area = numpy.pi * R**2
    u_l = G_l / area

    R_const = 8.314
    u_g = calculate_gas_velocity(
        G_gas=G_gas,
        p_t=p_t,
        T=T,
        R=R,
    )

    if u_g <= 0.0:
        raise ValueError("The calculated gas velocity must be positive.")

    c_in = K_S * numpy.sqrt(pl_in)
    c_g_in = pg_in / (R_const * T)

    # For LM, the isotope is atomic in the liquid and molecular in the gas.
    gas_to_liquid_ratio = 2.0 * u_g / u_l

    c_out_max_reaction = c_in - kla * Z / u_l * (c_in - K_S * numpy.sqrt(pg_in))

    # Maximum gas loading if the gas reaches equilibrium with the liquid inlet.
    c_g_max = pl_in / (R_const * T)

    c_out_max_gas = c_in - gas_to_liquid_ratio * (c_g_max - c_g_in)

    c_out_min = max(
        c_out_max_gas,
        c_out_max_reaction,
        K_S * numpy.sqrt(pg_in),
        0.0,
    )

    concentration_tolerance = 1.0e-12 * max(1.0, abs(c_in))

    if c_out_min > c_in + concentration_tolerance:
        raise ValueError(
            "No physical LM GLC extraction solution exists: "
            f"c_out_min={c_out_min:.6e} is greater than "
            f"c_in={c_in:.6e}."
        )

    def length_from_cout(c_out):
        liquid_outlet_pressure = c_out**2 / K_S**2

        return length_extractor_lm(
            R=R,
            G_l=G_l,
            G_gas=G_gas,
            pl_in=pl_in,
            pl_out=liquid_outlet_pressure,
            T=T,
            p_t=p_t,
            K_S=K_S,
            pg_in=pg_in,
            kla=kla,
            c_max=c_out_min,
        )

    # This is the largest height achievable before reaching the
    # gas-saturation/equilibrium/mass-transfer outlet limit.
    length_at_limit = length_from_cout(c_out_min)

    if length_at_limit <= Z + 1.0e-10:
        # The column is limited by the physical lower outlet bound.
        c_out = c_out_min
    else:

        def residual(c_out):
            return length_from_cout(c_out) - Z

        c_out = brentq(
            residual,
            c_out_min,
            c_in,
            xtol=1.0e-12,
            rtol=1.0e-12,
            maxiter=200,
        )

    c_out = float(numpy.clip(c_out, c_out_min, c_in))
    eff = 1.0 - c_out / c_in if c_in > 0.0 else 0.0

    return c_out, eff


def get_c_out_GLC_ms(
    Z,
    R,
    G_l,
    G_gas,
    pl_in,
    T,
    p_t,
    K_H,
    pg_in,
    kla,
):
    """
    Calculate the MS GLC outlet concentration and extraction efficiency.

    The molten-salt concentration is molecular Q2 concentration and
    therefore follows Henry's law, c_l = K_H * p_l.
    """

    if Z < 0.0:
        raise ValueError("The GLC height Z must be non-negative.")
    if R <= 0.0:
        raise ValueError("The GLC radius R must be positive.")
    if G_l <= 0.0:
        raise ValueError("The liquid flow rate G_l must be positive.")
    if G_gas <= 0.0:
        raise ValueError("The gas flow rate G_gas must be positive.")
    if kla < 0.0:
        raise ValueError("The mass-transfer coefficient kla cannot be negative.")

    area = numpy.pi * R**2
    u_l = G_l / area

    R_const = 8.314
    u_g = calculate_gas_velocity(
        G_gas=G_gas,
        p_t=p_t,
        T=T,
        R=R,
    )

    if u_g <= 0.0:
        raise ValueError("The calculated gas velocity must be positive.")

    c_in = K_H * pl_in
    c_g_in = pg_in / (R_const * T)

    # For MS, both liquid and gas concentrations refer to molecular Q2.
    gas_to_liquid_ratio = u_g / u_l

    c_out_max_reaction = c_in - kla * Z / u_l * (c_in - K_H * pg_in)

    # Maximum gas loading if the gas reaches equilibrium with the liquid inlet.
    c_g_max = pl_in / (R_const * T)

    c_out_max_gas = c_in - gas_to_liquid_ratio * (c_g_max - c_g_in)

    c_out_min = max(
        c_out_max_gas,
        c_out_max_reaction,
        K_H * pg_in,
        0.0,
    )

    concentration_tolerance = 1.0e-12 * max(1.0, abs(c_in))

    if c_out_min > c_in + concentration_tolerance:
        raise ValueError(
            "No physical MS GLC extraction solution exists: "
            f"c_out_min={c_out_min:.6e} is greater than "
            f"c_in={c_in:.6e}."
        )

    def length_from_cout(c_out):
        liquid_outlet_pressure = c_out / K_H

        # Important: this must call length_extractor_ms(), not
        # length_extractor_lm().
        return length_extractor_ms(
            R=R,
            G_l=G_l,
            G_gas=G_gas,
            pl_in=pl_in,
            pl_out=liquid_outlet_pressure,
            T=T,
            p_t=p_t,
            K_H=K_H,
            pg_in=pg_in,
            kla=kla,
            c_max=c_out_min,
        )

    # Largest achievable height before reaching the physical outlet limit.
    length_at_limit = length_from_cout(c_out_min)

    if length_at_limit <= Z + 1.0e-10:
        c_out = c_out_min
    else:

        def residual(c_out):
            return length_from_cout(c_out) - Z

        c_out = brentq(
            residual,
            c_out_min,
            c_in,
            xtol=1.0e-12,
            rtol=1.0e-12,
            maxiter=200,
        )
    c_out = float(numpy.clip(c_out, c_out_min, c_in))
    eff = 1.0 - c_out / c_in if c_in > 0.0 else 0.0

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
