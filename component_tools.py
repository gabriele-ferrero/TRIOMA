import numpy as np
import molten_salts as MS
import liquid_metals as LM
import correlations as corr
import extractor
import matplotlib.pyplot as plt


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
        eff (float): The efficiency of the component.
        fluid (Fluid, optional): The fluid associated with the component. Defaults to None.
        membrane (Membrane, optional): The membrane associated with the component. Defaults to None.
    """

    def __init__(
        self,
        c_in: float,
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

    def outlet_c_comp(self) -> float:
        """
        Calculates the concentration of the component at the outlet.

        Returns:
            float: The concentration of the component at the outlet.
        """
        c_out = self.c_in * (1 - self.eff)
        return c_out

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
                        P_H2=self.fluid.p_H2,
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
                        P_H2=self.fluid.p_H2,
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
                    P_H2=self.fluid.p_H2,
                )
            else:
                self.H = LM.partition_param(
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
                    P_H2=self.fluid.p_H2,
                )

    def get_efficiency(self, L):  # L is the characteristic length of the component
        """
        Calculates the efficiency of the component.

        Args:
            L (float): The characteristic length of the component.

        Returns:
            float: The efficiency of the component.
        """
        L_vec = np.linspace(0, L, 1000)
        dl = L_vec[1] - L_vec[0]

        c_vec = np.ndarray(len(L_vec))
        for i in range(len(L_vec)):
            if self.fluid.MS:
                f_H2 = 0.5
            else:
                f_H2 = 1
            if i == 0:

                c_vec[i] = self.c_in
                self.get_flux(c_vec[i])
            else:
                c_vec[i] = c_vec[
                    i - 1
                ] + f_H2 * self.J_perm * self.fluid.d_Hyd * np.pi * dl**2 / self.fluid.U0 / (
                    np.pi * self.fluid.d_Hyd**2 / 4 * dl
                )
                self.get_flux(c_vec[i])
        plt.plot(L_vec, c_vec)
        self.eff = 1 - (self.c_in - c_vec[-1]) / self.c_in

    def get_flux(self, c):
        self.get_adimensionals()
        if self.fluid.MS:
            if self.H < 0.1:  ## TODO limits to check
                # print("Surface limited approximation")
                self.J_perm = -self.membrane.k_d * (c / self.fluid.Solubility)
            elif self.W > 10:  ## TODO limits to check
                # print("Diffusion limited approximation")
                self.J_perm = -(
                    self.membrane.D
                    / self.membrane.thick
                    * self.membrane.K_S
                    * (c / self.fluid.Solubility) ** 0.5
                )
            elif self.H > 10 and self.H / self.W > 100:  ## TODO limits to check
                # print("Mass transport limited approximation")

                self.J_perm = -2 * self.fluid.k_t * c
            else:
                print("Mixed regime approximation")
                raise ValueError("Mixed regime not yet implemented")
        else:  ### Liquid Metal
            raise ValueError("Liquid Metal not implemented")
            if self.H < 0.1:  ## TODO limits to check
                # print("Surface limited approximation")
                self.J_perm = -self.membrane.k_r * (c**2)
            elif self.W > 10:  ## TODO limits to check
                # print("Diffusion limited approximation")
                self.J_perm = -(
                    self.membrane.D
                    / self.membrane.thick
                    * self.membrane.K_S
                    * (c / self.fluid.Solubility)
                )
            elif self.H > 10 and self.H / self.W > 100:  ## TODO limits to check
                print("Mass transport limited approximation")

                self.J_perm = -self.fluid.k_t * c
            else:
                print("Mixed regime approximation")
                raise ValueError("Mixed regime not implemented")


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
        c0 (float, optional): Inlet Tritium Concentration of the fluid. Defaults to None.
        p_H2 (float, optional): Hydrogen pressure of the fluid at the inlet. Defaults to None.
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
        c0: float,
        d_Hyd: float = None,
        k_t: float = None,
        mu: float = None,
        rho: float = None,
        U0: float = None,
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
            c0 (float, optional): Inlet Concentration of the fluid. Defaults to None.
            p_H2 (float, optional): Hydrogen pressure of the fluid at the inlet. Defaults to None.
            mu (float, optional): Viscosity of the fluid. Defaults to None.
            rho (float, optional): Density of the fluid. Defaults to None.
            U0 (float, optional): Velocity of the fluid. Defaults to None.
        """
        self.T = T
        self.Solubility = Solubility
        self.MS = MS
        self.D = D
        self.k_t = k_t
        self.c0 = c0
        if self.MS:
            self.p_H2 = c0 / Solubility
        else:
            self.p_H2 = (c0 / Solubility) ** 2
        self.d_Hyd = d_Hyd
        self.mu = mu
        self.rho = rho
        self.U0 = U0

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

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
                self.k_t = corr.get_k_from_Sh(
                    corr.Sherwood(
                        corr.Schmidt(self.D, self.mu, self.rho),
                        corr.Re(self.rho, self.U0, self.d_Hyd, self.mu),
                    ),
                    self.d_Hyd,
                    self.D,
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

    def update_attribute(self, attr_name, new_value):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)


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
