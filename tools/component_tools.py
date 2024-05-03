import numpy as np

# from example_simulation import TBR
import tools.molten_salts as MS
import tools.liquid_metals as LM
import tools.correlations as corr
import matplotlib.pyplot as plt
import tools.extractor as extractor
from scipy.constants import N_A
from scipy.constants import physical_constants
from scipy.optimize import minimize
from scipy.special import lambertw


def print_class_variables(instance):
    """
    Prints all variables of a class. If a variable is a class itself,
    calls this function recursively for the internal class.

    Args:
        instance: The class instance.
    """
    built_in_types = [
        Component,
        Fluid,
        Membrane,
        GLC_Gas,
        GLC,
        FluidMaterial,
        SolidMaterial,
        BreedingBlanket,
    ]
    for attr_name, attr_value in instance.__dict__.items():
        if type(attr_value) in built_in_types:
            print(f"{attr_name} is a {type(attr_value)} class, printing its variables:")
            print_class_variables(attr_value)
        else:
            print(f"{attr_name}: {attr_value}")


def calculate_p_H2_from_c0(instance, c0):
    return c0 / instance.fluid.Solubility


def set_attribute(instance, attr_name, new_value):
    """
    Sets the specified attribute to a new value.

    Args:
        instance: The class instance.
        attr_name (str): The name of the attribute to set.
        new_value: The new value for the attribute.
    """
    if hasattr(instance, attr_name):
        setattr(instance, attr_name, new_value)
    else:
        raise ValueError(
            f"'{attr_name}' is not an attribute of {instance.__class__.__name__}"
        )


class Component:
    """
    Represents a component in a plant to make a high level T transport analysis.

    Args:
        c_in (float): The concentration of the component at the inlet.
        fluid (Fluid, optional): The fluid associated with the component. Defaults to None.
        membrane (Membrane, optional): The membrane associated with the component. Defaults to None.
    """

    def __init__(
        self,
        c_in: float = None,
        eff: float = None,
        fluid: "Fluid" = None,
        membrane: "Membrane" = None,
    ):
        """
        Initializes a new instance of the Component class.

        Args:
            c_in (float): The concentration of the component at the inlet.
            eff (float, optional): The efficiency of the component. Defaults to None.
            fluid (Fluid, optional): The fluid associated with the component. Defaults to None.
            membrane (Membrane, optional): The membrane associated with the component. Defaults to None.
        """
        self.c_in = c_in
        self.eff = eff
        self.fluid = fluid
        self.membrane = membrane
        self.H = None
        self.W = None
        if self.fluid.k_t is None:
            self.fluid.get_kt()

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)

    def outlet_c_comp(self) -> float:
        """
        Calculates the concentration of the component at the outlet.

        Returns:
            float: The concentration of the component at the outlet.
        """
        self.c_out = self.c_in * (1 - self.eff)

    def T_leak(self) -> float:
        """
        Calculates the leakage of the component.

        Returns:
            float: The leakage of the component.
        """
        leak = self.c_in * self.eff
        return leak

    def get_regime(self):
        """
        Gets the regime of the component.

        Returns:
            str: The regime of the component.
        """
        if self.fluid is not None:

            if self.fluid.MS:
                if self.membrane is not None:
                    return MS.get_regime(
                        k_d=self.membrane.k_d,
                        D=self.membrane.D,
                        thick=self.membrane.thick,
                        K_S=self.membrane.K_S,
                        c0=self.c_in,
                        k_t=self.fluid.k_t,
                        k_H=self.fluid.Solubility,
                    )
                else:
                    return "No membrane selected"
            else:
                if self.membrane is not None:
                    return LM.get_regime(
                        D=self.membrane.D,
                        k_t=self.fluid.k_t,
                        K_S_S=self.membrane.K_S,
                        K_S_L=self.fluid.Solubility,
                        k_r=self.membrane.k_r,
                        thick=self.membrane.thick,
                        c0=self.c_in,
                    )
                else:
                    return "No membrane selected"

    def get_adimensionals(self):
        """
        Calculates the adimensional parameters H and W.

        Updates the H and W attributes of the Component object.
        """
        if self.fluid is not None:
            if self.fluid.MS:
                self.H = MS.H(
                    k_t=self.fluid.k_t, k_H=self.fluid.Solubility, k_d=self.membrane.k_d
                )
                self.W = MS.W(
                    k_d=self.membrane.k_d,
                    D=self.membrane.D,
                    thick=self.membrane.thick,
                    K_S=self.membrane.K_S,
                    c0=self.c_in,
                    k_H=self.fluid.Solubility,
                )
            else:
                self.H = LM.W(
                    k_r=self.membrane.k_r,
                    D=self.membrane.D,
                    thick=self.membrane.thick,
                    K_S=self.membrane.K_S,
                    c0=self.c_in,
                    K_S_L=self.fluid.Solubility,
                ) * LM.partition_param(
                    D=self.membrane.D,
                    k_t=self.fluid.k_t,
                    K_S_S=self.membrane.K_S,
                    K_S_L=self.fluid.Solubility,
                    t=self.membrane.thick,
                )
                self.W = LM.W(
                    k_r=self.membrane.k_r,
                    D=self.membrane.D,
                    thick=self.membrane.thick,
                    K_S=self.membrane.K_S,
                    c0=self.c_in,
                    K_S_L=self.fluid.Solubility,
                )

    def use_analytical_efficiency(self, L):
        """Evaluates the analytical efficiency and substitutes it in the efficiency attribute of the component.

        Args:
            L (float): the length of the pipe component
        Returns:
            None
        """
        self.analytical_efficiency(L)
        self.eff = self.eff_an

    def get_efficiency(self, L, plotvar: bool = False, c_guess: float = None):
        """
        Calculates the efficiency of the component.

        Args:
            L (float): The characteristic length of the component.

        Returns:
            float: The efficiency of the component.
        """
        L_vec = np.linspace(0, L, 100)
        dl = L_vec[1] - L_vec[0]

        c_vec = np.ndarray(len(L_vec))
        for i in range(len(L_vec)):
            if self.fluid.MS:
                f_H2 = 0.5
            else:
                f_H2 = 1
            if i == 0:

                c_vec[i] = self.c_in

                if c_guess is float:
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=None)
            else:
                c_vec[i] = c_vec[
                    i - 1
                ] + f_H2 * self.J_perm * self.fluid.d_Hyd * np.pi * dl**2 / self.fluid.U0 / (
                    np.pi * self.fluid.d_Hyd**2 / 4 * dl
                )
                if c_guess is float:
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=None)
        if plotvar:
            plt.plot(L_vec, c_vec)
        self.eff = (self.c_in - c_vec[-1]) / self.c_in

    def analytical_efficiency(self, L):
        """
        Calculate the analytical efficiency of a component.

        Parameters:
        - L: Length of the component

        Returns:
        - eff_an: Analytical efficiency of the component (from Humrickhouse papers)

        This function calculates the analytical efficiency of a component based on the given length (L) of the component.
        It uses various properties of the component, such as membrane properties, fluid properties, and adimensionals.



        If the fluid is a molten salt (MS=True), the analytical efficiency is calculated using the following formula:
        eff_an = 1 - epsilon * (lambertw(z=np.exp(beta - tau - 1), tol=1e-10) ** 2 + 2 * lambertw(z=np.exp(beta - tau - 1), tol=1e-10)).

        If the fluid is a liquid metal (MS= False), the analytical efficiency is calculated using the following formula:
        eff_an = 1 - np.exp(-tau * zeta / (1 + zeta)).

        The output of the function is the analytical efficiency of the component as Component.eff_an.
        """

        self.get_adimensionals()
        self.tau = 4 * self.fluid.k_t * L / (self.fluid.U0 * self.fluid.d_Hyd)
        if self.fluid.MS:
            epsilon = (
                1
                / self.c_in
                / self.fluid.Solubility
                * (
                    0.5
                    * self.membrane.K_S
                    * self.membrane.D
                    / (
                        self.fluid.k_t
                        * self.fluid.d_Hyd
                        * np.log(
                            (self.fluid.d_Hyd + 2 * self.membrane.thick)
                            / self.fluid.d_Hyd
                        )
                    )
                )
                ** 2
            )
            beta = (1 / epsilon + 1) ** 0.5 + np.log((1 / epsilon + 1) ** 0.5 - 1)
            max_exp = np.log(np.finfo(np.float64).max)
            beta_tau = beta - self.tau - 1
            if beta_tau > max_exp:
                print("Warning: Overflow encountered in exp, input too large.")
                # Handle the overflow case here, e.g., by setting a maximum value
                z = np.finfo(np.float64).max
            else:
                z = np.exp(beta_tau)
            w = lambertw(z, tol=1e-10)
            self.eff_an = 1 - epsilon * (w**2 + 2 * w)
            if self.eff_an.imag != 0:
                raise ValueError("self.eff_an has a non-zero imaginary part")
            else:
                self.eff_an = self.eff_an.real  # get rid of 0*j
        else:
            self.zeta = (2 * self.membrane.K_S * self.membrane.D) / (
                self.fluid.k_t
                * self.fluid.Solubility
                * self.fluid.d_Hyd
                * np.log(
                    (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                )
            )
            self.eff_an = 1 - np.exp(-self.tau * self.zeta / (1 + self.zeta))

    def get_flux(self, c, c_guess=None):
        """
        Calculates the Tritium flux of the component.
        It can make some approximations based on W and H to make the solver faster

        Args:
            c (float): The concentration in the component fluid.

        Returns:
            float: The permeation flux.

        """
        self.get_adimensionals()
        if self.fluid.MS:
            if self.W > 10:
                # DIFFUSION LIMITED
                if self.H / self.W > 1000:
                    # Mass transport limited
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
                elif self.H / self.W < 0.0001:
                    # Diffusion limited
                    self.J_perm = -(
                        self.membrane.D
                        / (
                            self.fluid.d_Hyd
                            / 2
                            * np.log(
                                (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                / (self.fluid.d_Hyd / 2)
                            )
                        )
                        * self.membrane.K_S
                        * (c / self.fluid.Solubility) ** 0.5
                    )
                else:
                    # Mixed regime mass transport diffusion
                    def equations(vars):
                        c_wl = vars
                        J_mt = 2 * self.fluid.k_t * (c - c_wl)
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (
                                self.membrane.K_S
                                * (c_wl / self.fluid.Solubility) ** 0.5
                            )
                        )
                        return abs((J_diff - J_mt))

                    if c_guess is None:
                        initial_guess = [(c * 1e-2)]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    self.J_perm = (
                        -2 * self.fluid.k_t * (c - solution.x[0])
                    )  ## MS factor
                    return float(solution.x[0])
            elif self.W < 0.1:
                # Surface limited
                if self.H > 100:
                    # Mass transport limited
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
                elif self.H < 0.01:
                    # Surface limited
                    self.J_perm = -self.membrane.k_d * (c / self.fluid.Solubility)
                else:
                    # Mixed regime mass transfer surface
                    def equations(vars):
                        c_wl = vars
                        c_bl = c
                        J_mt = 2 * self.fluid.k_t * (c_bl - c_wl)  ## MS factor
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )

                        return abs(J_mt - J_surf)

                    initial_guess = [(c * 1e-1)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_bl = c
                    c_wl = solution.x[0]
                    self.J_perm = 2 * self.fluid.k_t * (c_bl - c_wl)  ## MS factor
                    return float(solution.x[0])
            else:
                if self.H / self.W < 0.0001:
                    # Mass transport limited
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
                elif self.H / self.W > 1000:
                    # Mixed Diffusion Surface
                    def equations(vars):
                        c_wl = vars
                        c_bl = c
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (
                                self.membrane.K_S
                                * (c_wl / self.fluid.Solubility) ** 0.5
                            )
                        )

                        return abs(J_diff - J_surf)

                    if c_guess is None:
                        initial_guess = [(c / 2)]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (1e-14, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_wall = solution.x[0]
                    self.J_perm = (
                        self.membrane.D
                        / (
                            self.fluid.d_Hyd
                            / 2
                            * np.log(
                                (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                / (self.fluid.d_Hyd / 2)
                            )
                        )
                        * (self.membrane.K_S * (c_wall / self.fluid.Solubility) ** 0.5)
                    )
                    return float(solution.x[0])
                else:
                    # Mixed regime mass transport diffusion surface and diffusion
                    def equations(vars):
                        c_wl, c_ws = vars

                        c_bl = c
                        J_mt = 2 * self.fluid.k_t * (c_bl - c_wl)  ## MS factor

                        J_d = self.membrane.k_d * (
                            c_wl / self.membrane.K_S
                        ) - self.membrane.k_d * self.membrane.K_S * (c_ws**2)
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (self.membrane.K_S * c_ws)
                        )
                        eq1 = abs(J_mt / J_d - 1)
                        eq2 = abs(J_mt / J_diff - 1)

                        return [eq1, eq2]

                    initial_guess = [(2 * c / 3), (c / 3)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_wl = solution.x[0]
                    self.J_perm = 2 * self.fluid.k_t * (c - c_wl)
                    return float(solution.x[0])
        else:
            if self.W > 10:
                # DIFFUSION LIMITED
                if self.H / self.W > 1000:
                    # Mass transport limited
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
                elif self.H / self.W < 0.0001:
                    # Diffusion limited
                    self.J_perm = -(
                        self.membrane.D
                        / (
                            self.fluid.d_Hyd
                            / 2
                            * np.log(
                                (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                / (self.fluid.d_Hyd / 2)
                            )
                        )
                        * self.membrane.K_S
                        * (c / self.fluid.Solubility)
                    )
                else:
                    # Mixed regime mass transport diffusion
                    def equations(vars):
                        c_wl = vars
                        J_mt = self.fluid.k_t * (c - c_wl)
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (self.membrane.K_S * (c_wl / self.fluid.Solubility))
                        )
                        return abs((J_diff - J_mt))

                    if c_guess is None:
                        initial_guess = [(c * 1e-2)]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    self.J_perm = -self.fluid.k_t * (c - solution.x[0])  ## LM factor
                    return float(solution.x[0])
            elif self.W < 0.1:
                # Surface limited
                if self.H > 100:
                    # Mass transport limited
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
                elif self.H < 0.01:
                    # Surface limited
                    self.J_perm = -self.membrane.k_d * (c / self.fluid.Solubility)
                else:
                    # Mixed regime mass transfer surface
                    def equations(vars):
                        c_wl = vars
                        c_bl = c
                        J_mt = self.fluid.k_t * (c_bl - c_wl)  ## LM factor
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )

                        return abs(J_mt - J_surf)

                    initial_guess = [(c * 1e-1)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_bl = c
                    c_wl = solution.x[0]
                    self.J_perm = self.fluid.k_t * (c_bl - c_wl)  ## LM factor
                    return float(solution.x[0])
            else:
                if self.H / self.W < 0.0001:
                    # Mass transport limited
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
                elif self.H / self.W > 1000:
                    # Mixed Diffusion Surface
                    def equations(vars):
                        c_wl = vars
                        c_bl = c
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (self.membrane.K_S * (c_wl / self.fluid.Solubility))
                        )

                        return abs(J_diff - J_surf)

                    if c_guess is None:
                        initial_guess = [(c / 2)]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (1e-14, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_wall = solution.x[0]
                    self.J_perm = (
                        self.membrane.D
                        / (
                            self.fluid.d_Hyd
                            / 2
                            * np.log(
                                (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                / (self.fluid.d_Hyd / 2)
                            )
                        )
                        * (self.membrane.K_S * (c_wall / self.fluid.Solubility))
                    )
                    return float(solution.x[0])
                else:
                    # Mixed regime mass transport diffusion surface and diffusion
                    def equations(vars):
                        c_wl, c_ws = vars

                        c_bl = c
                        J_mt = self.fluid.k_t * (c_bl - c_wl)  ## LM factor

                        J_d = self.membrane.k_d * (
                            c_wl / self.membrane.K_S
                        ) - self.membrane.k_d * self.membrane.K_S * (c_ws**2)
                        J_diff = (
                            self.membrane.D
                            / (
                                self.fluid.d_Hyd
                                / 2
                                * np.log(
                                    (self.fluid.d_Hyd / 2 + self.membrane.thick)
                                    / (self.fluid.d_Hyd / 2)
                                )
                            )
                            * (self.membrane.K_S * c_ws)
                        )
                        eq1 = abs(J_mt / J_d - 1)
                        eq2 = abs(J_mt / J_diff - 1)

                        return [eq1, eq2]

                    initial_guess = [(2 * c / 3), (c / 3)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c),
                        ],
                        tol=1e-15,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_wl = solution.x[0]
                    self.J_perm = self.fluid.k_t * (c - c_wl)  # LM factor
                    return float(solution.x[0])

    def get_global_HX_coeff(self, R_conv_sec: float = 0):
        """
        Calculates the global heat exchange coefficient of the component.
        It can take the secondary resistance to convection as input. defaults to no resistance

        Returns:
            float: The global heat exchange coefficient of the component.
        """
        R_cond = np.log((self.fluid.d_Hyd + self.membrane.thick) / self.fluid.d_Hyd) / (
            2 * np.pi * self.membrane.k
        )
        Re = corr.Re(self.fluid.rho, self.fluid.U0, self.fluid.d_Hyd, self.fluid.mu)
        Pr = corr.Pr(self.fluid.cp, self.fluid.mu, self.fluid.k)
        h_prim = corr.get_h_from_Nu(
            corr.Nu_DittusBoelter(Re, Pr), self.fluid.k, self.fluid.d_Hyd
        )
        R_conv_prim = 1 / h_prim
        R_tot = R_conv_prim + R_cond + R_conv_sec
        self.U = 1 / R_tot
        return


class Fluid:
    """
    Represents a fluid in a component for Tritium transport analysis

    Args:
        T (float): Temperature of the fluid.
        D (float): Tritium Diffusivity of the fluid.
        Solubility (float): Solubility of the fluid.
        MS (bool): Indicates whether the fluid is a molten salt or a liquid metal.
        d_Hyd (float, optional): Hydraulic diameter of the fluid. Defaults to None.
        k_t (float, optional): Mass transport coefficient of the fluid. Defaults to None.
        mu (float, optional): Viscosity of the fluid. Defaults to None.
        rho (float, optional): Density of the fluid. Defaults to None.
        U0 (float, optional): Velocity of the fluid. Defaults to None.
    """

    def __init__(
        self,
        T: float,
        D: float,
        Solubility: float,
        MS: bool,
        d_Hyd: float = None,
        k_t: float = None,
        mu: float = None,
        rho: float = None,
        U0: float = None,
        k: float = None,
        cp: float = None,
    ):
        """
        Initializes a new instance of the Fluid class.

        Args:
            T (float): Temperature of the fluid.
            D (float): T Diffusivity of the fluid.
            Solubility (float): Solubility of the fluid.
            MS (bool): Indicates whether the fluid is a molten salt or a liquid metal.
            d_Hyd (float, optional): Hydraulic diameter of the fluid. Defaults to None.
            k_t (float, optional): Mass transport coefficient of the fluid. Defaults to None.
            mu (float, optional): Viscosity of the fluid. Defaults to None.
            rho (float, optional): Density of the fluid. Defaults to None.
            U0 (float, optional): Velocity of the fluid. Defaults to None.
            k thermal conductivity of the fluid. Defaults to None.
        """
        self.T = T
        self.Solubility = Solubility
        self.MS = MS
        self.D = D
        self.k_t = k_t
        self.d_Hyd = d_Hyd
        self.mu = mu
        self.rho = rho
        self.U0 = U0
        self.k = k
        self.cp = cp

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)

    def get_kt(self):
        """
        Calculates the mass transport coefficient (k_t) for the fluid.

        If the hydraulic diameter (d_Hyd) is defined, the mass transport coefficient is calculated using correlations.
        Otherwise, an error message is printed.

        Returns:
            None
        """
        if self.d_Hyd:
            if self.k_t is None:
                Re = corr.Re(rho=self.rho, u=self.U0, L=self.d_Hyd, mu=self.mu)
                Sc = corr.Schmidt(D=self.D, mu=self.mu, rho=self.rho)
                if Re < 1e4 and Re > 2030:
                    Sh = 0.015 * Re**0.83 * Sc**0.42  ## Stempien Thesis pg 155-157
                elif Re > 1e4:
                    Sh = 0.0096 * Re**0.913 * Sc**0.346  ##Getthem paper
                else:
                    print(Re)
                    raise ValueError("Reynolds number is too low")
                self.k_t = corr.get_k_from_Sh(
                    Sh=Sh,
                    L=self.d_Hyd,
                    D=self.D,
                )

            else:
                print("k_t is already defined")
        else:
            print("Hydraulic Diameter is not defined")


class Membrane:
    """
    Represents a metallic membrane of a component for H transport.

    Attributes:
        T (float): Temperature of the membrane.
        D (float): Diffusion coefficient of the membrane.
        thick (float): Thickness of the membrane.
        K_S (float): Solubility coefficient of the membrane.
        k_d (float, optional): Dissociation rate constant of the membrane. Defaults to 1e6.
        k_r (float, optional): Recombination rate constant of the membrane. Defaults to 1e6.
    """

    def __init__(
        self,
        T: float,
        D: float,
        thick: float,
        K_S: float,
        k_d: float = 1e6,
        k_r: float = 1e6,
        k: float = None,
    ):
        """
        Initializes a new instance of the Membrane class.

        Args:
            T (float): Temperature of the membrane.
            D (float): Diffusion coefficient of the membrane.
            thick (float): Thickness of the membrane.
            K_S (float): Solubility coefficient of the membrane.
            k_d (float, optional): Dissociation rate constant of the membrane. Defaults to 1e6.
            k_r (float, optional): Recombination rate constant of the membrane. Defaults to 1e6.
        """
        self.T = T
        self.D = D
        self.thick = thick
        self.k_d = k_d
        self.K_S = K_S
        self.k_r = k_r
        self.k = k

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)


class GLC_Gas:
    """
    GLC_Gas class represents a sweep gas in a GLC (Gas-Liquid Contactor) system.

    Attributes:
        G_gas (float): The flow rate of the gas.
        pg_in (float, optional): The Tritium inlet partial pressure of the gas. Default is 0.
        p_tot (float, optional): The total pressure of the component. Default is 100000 Pa.
        kla (float, optional): The total (kl*Area) mass transfer coefficient. Default is 0.
    """

    def __init__(
        self, G_gas: float, pg_in: float = 0, p_tot: float = 100000, kla: float = 0
    ):
        """
        Initializes a new instance of the GLC_Gas class.

        Args:
            G_gas (float): The flow rate of the gas.
            pg_in (float, optional): The Tritium inlet partial pressure of the gas. Default is 0.
            p_tot (float, optional): The total pressure of the component. Default is 100000 Pa.
            kla (float, optional): The total (kl*Area) mass transfer coefficient. Default is 0.
        """
        self.G_gas = G_gas
        self.pg_in = pg_in
        self.p_tot = p_tot
        self.kla = kla

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)


class GLC(Component):
    """
    GLC (Gas-Liquid Contact) class represents a gas-liquid contactor component.

    Args:
        H (float): Height of the GLC.
        R (float): Radius of the GLC.
        L (float): Characteristic Length of the GLC fluid flow.
        c_in (float): Inlet concentration of the GLC.
        eff (float, optional): Efficiency of the GLC. Defaults to None.
        fluid (Fluid, optional): Fluid object representing the liquid phase. Defaults to None.
        membrane (Membrane, optional): Membrane object representing the membrane used in the GLC. Not super important. Defaults to None.
        GLC_gas (GLC_Gas, optional): GLC_Gas object representing the gas phase. Defaults to None.

    Attributes:
        H (float): Height of the GLC [m].
        R (float): Radius of the GLC [m].
        L (float): Length of the GLC [m].
        GLC_gas (GLC_Gas): GLC_Gas object representing the gas phase.

    Methods:
        get_kla_Ring(): Calculates the mass transfer coefficient (kla) for a Raschig Ring matrix.
    """

    def __init__(
        self,
        H: float,
        R: float,
        L: float,
        c_in: float,
        eff: float = None,
        fluid: "Fluid" = None,
        membrane: "Membrane" = None,
        GLC_gas: "GLC_Gas" = None,
    ):
        """
        Initializes a new instance of the GLC class.

        Args:
            H (float): Height of the GLC.
            R (float): Radius of the GLC.
            L (float): Characteristic Length of the GLC fluid flow.
            c_in (float): Inlet concentration of the GLC.
            eff (float, optional): Efficiency of the GLC. Defaults to None.
            fluid (Fluid, optional): Fluid object representing the liquid phase. Defaults to None.
            membrane (Membrane, optional): Membrane object representing the membrane used in the GLC. Not super important. Defaults to None.
            GLC_gas (GLC_Gas, optional): GLC_Gas object representing the gas phase. Defaults to None.
        """
        super().__init__(c_in, eff, fluid, membrane)
        self.H = H
        self.R = R
        self.L = L
        self.GLC_gas = GLC_gas

    def get_kla_Ring(self):
        """
        Calculates the mass transfer coefficient (kla) for a Raschig ring matrix.

        The mass transfer coefficient is calculated based on the Reynolds number (Re) and Schmidt number (Sc) of the fluid,
        as well as the diameter (d) of the ring.

        Returns:
            None
        """
        d = 2e-3  # Ring diameter
        Re = corr.Re(rho=self.fluid.rho, u=self.fluid.U0, L=self.L, mu=self.fluid.mu)
        Sc = corr.Schmidt(D=self.fluid.D, mu=self.fluid.mu, rho=self.fluid.rho)
        self.GLC_gas.kla = extractor.corr_packed(
            Re,
            Sc,
            d,
            rho_L=self.fluid.rho,
            mu_L=self.fluid.mu,
            L=self.L,
            D=self.fluid.D,
        )

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)


class FluidMaterial:
    """
    Represents a fluid material with various properties.

    Attributes:
        T (float): Temperature of the fluid material.
        D (float): Density of the fluid material.
        Solubility (float): Solubility of the fluid material.
        MS (float): Molecular weight of the fluid material.
        mu (float): Viscosity of the fluid material.
        rho (float): Density of the fluid material.
        k (float): Thermal conductivity of the fluid material.
        cp (float): Specific heat capacity of the fluid material.
    """

    def __init__(self, T, D, Solubility, MS, mu, rho, k, cp):
        self.T = T
        self.D = D
        self.Solubility = Solubility
        self.MS = MS
        self.mu = mu
        self.rho = rho
        self.k = k
        self.cp = cp

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)


class SolidMaterial:
    """
    Represents a solid material used in a component.

    Attributes:
        D (float): The Diffusivity of the solid material.
        K_S (float): The Sievert constant of the solid material.
        k_d (float): Dissociation constant of the solid material.
        k_r (float): The Recombination constant of the solid material.
    """

    def __init__(self, D, K_S, k_d, k_r):
        self.D = D
        self.K_S = K_S
        self.k_d = k_d
        self.k_r = k_r

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)


class BreedingBlanket:
    """
    Represents a breeding blanket component in a fuel cycle system.

    Attributes:
        c_in(float): Inlet concentration of tritium in the breeding blanket component.
        Q (float): Heat generated by the breeding blanket component.
        TBR (float): Tritium breeding ratio of the breeding blanket component.
        T_out (float): Outlet temperature of the breeding blanket component.
        T_in (float): Inlet temperature of the breeding blanket component.
        fluid (Fluid): Fluid used in the breeding blanket component.
    """

    def __init__(
        self,
        c_in: float = None,
        Q: float = None,
        TBR: float = None,
        T_out: float = None,
        T_in: float = None,
        fluid: Fluid = None,
    ):
        self.c_in = c_in
        self.Q = Q
        self.TBR = TBR
        self.T_out = T_out
        self.T_in = T_in
        self.fluid = fluid

    def inspect(self):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self)

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def get_flowrate(self):
        """
        Calculates the flow rate of the coolant in the breeding blanket component.
        """
        self.m_coolant = self.Q / ((self.T_out - self.T_in) * self.fluid.cp)
        return

    def get_cout(self, print_var: bool = False):
        """
        Calculates the outlet concentration of tritium in the breeding blanket component.

        Args:
            print_var (bool): If True, prints the intermediate variables.

        Returns:
            None
        """
        self.get_flowrate()
        eV_to_J = physical_constants["electron volt-joule relationship"][0]
        reaction_energy = 17.6e6  # reaction energy in eV 17.6 MeV
        neutrons = self.Q / (reaction_energy * eV_to_J)

        tritium_gen = self.TBR * neutrons / N_A  ##moles
        if print_var:
            print("neu", neutrons)
            print("Trit", tritium_gen)
        self.c_out = tritium_gen / (self.m_coolant / self.fluid.rho) + self.c_in
