import numpy as np
from typing import Union

from TRIOMA.tools import correlations as corr
from TRIOMA.tools.TriomaClass import TriomaClass


class Geometry(TriomaClass):
    """
    Represents the geometry of a component.

    Args:
        L (float): Length of the component.
        D (float): Diameter of the component.
        thick (float): Thickness of the component.
        n_pipes (int, optional): The number of pipes in the component. Defaults to 1.

    """

    def __init__(
        self,
        L: float = None,
        D: float = None,
        thick: float = None,
        n_pipes: float = 1,
        turbulator: Union["Turbulator"] = None,
    ):
        self.L = L
        self.D = D
        self.thick = thick
        self.n_pipes = n_pipes
        self.turbulator = turbulator

    def get_fluid_volume(self):
        """
        Calculates the volume of the fluid component.
        """

        return np.pi * (self.D / 2) ** 2 * self.L

    def get_solid_volume(self):
        """
        Calculates the volume of the solid component.
        """
        return np.pi * ((self.D / 2) ** 2 - (self.D / 2 - self.thick) ** 2) * self.L

    def get_total_volume(self):
        """
        Calculates the total volume of the component.
        """
        return self.get_fluid_volume() + self.get_solid_volume()


class Fluid(TriomaClass):
    """
    Represents a fluid in a component for Tritium transport analysis

    Args:
        T (float): Temperature of the fluid.
        D (float): Tritium Diffusivity of the fluid.
        D_0 (float): Preexponential Diffusivity of the fluid.
        E_d (float): Activation energy of the fluid diffusivity.
        Solubility (float): Solubility of the fluid.
        Solubility_0 (float): Preexponential solubility of the fluid.
        E_s (float): Activation energy of the fluid solubility.
        MS (bool): Indicates whether the fluid is a molten salt or a liquid metal.
        d_Hyd (float, optional): Hydraulic diameter of the fluid. Defaults to None.
        k_t (float, optional): Mass transport coefficient of the fluid. Defaults to None.
        mu (float, optional): Viscosity of the fluid. Defaults to None.
        rho (float, optional): Density of the fluid. Defaults to None.
        U0 (float, optional): Velocity of the fluid. Defaults to None.
        inv (float, optional): Inventory of the fluid. Defaults to None.
        recirculation(float,optional): fraction of recirculated flowrate. Defaults to 0. 1 = 100%, 0.5=50%.
    """

    def __init__(
        self,
        T: float = None,
        D: float = None,
        D_0: float = None,
        E_d: float = None,
        Solubility: float = None,
        Solubility_0: float = None,
        E_s: float = None,
        MS: bool = True,
        d_Hyd: float = None,
        k_t: float = None,
        mu: float = None,
        rho: float = None,
        U0: float = None,
        k: float = None,
        cp: float = None,
        inv: float = None,
        recirculation: float = 0,
        V: float = None,
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
        self.MS = MS
        self.D_0 = D_0
        self.E_d = E_d
        self.Solubility_0 = Solubility_0
        self.E_s = E_s
        if D_0 is not None and E_d is not None:
            self.D = D_0 * np.exp(-E_d / (8.617333262145e-5 * self.T))
        else:
            self.D = D
        if Solubility_0 is not None and E_s is not None:
            self.Solubility = Solubility_0 * np.exp(-E_s / (8.617333262145e-5 * self.T))
        else:
            self.Solubility = Solubility
        self.k_t = k_t
        self.d_Hyd = d_Hyd
        self.mu = mu
        self.rho = rho
        self.recirculation = recirculation
        if self.recirculation == 0:
            self.U0 = U0
        else:
            self.U0 = U0 * (1 + self.recirculation)
        self.k = k
        self.cp = cp
        self.inv = inv

    def set_properties_from_fluid_material(
        self, fluid_material: "FluidMaterial" = None
    ):
        """
        Sets the properties of the fluid from a FluidMaterial object.

        Args:
            fluid_material (FluidMaterial): The FluidMaterial object to set the properties from.
        """
        self.T = fluid_material.T
        self.D = fluid_material.D
        self.Solubility = fluid_material.Solubility
        self.mu = fluid_material.mu
        self.rho = fluid_material.rho
        self.cp = fluid_material.cp
        self.k = fluid_material.k

    def update_T_prop(self):
        if self.D_0 is not None and self.E_d is not None:
            self.D = self.D_0 * np.exp(-self.E_d / (8.617333262145e-5 * self.T))
        if self.Solubility_0 is not None and self.E_s is not None:
            self.Solubility = self.Solubility_0 * np.exp(
                -self.E_s / (8.617333262145e-5 * self.T)
            )
        return

    def get_kt(self, turbulator=None):
        """
        Calculates the mass transport coefficient (k_t) for the fluid.

        If the hydraulic diameter (d_Hyd) is defined, the mass transport coefficient is calculated using correlations.
        Otherwise, an error message is printed.

        Returns:
            None
        """
        if self.d_Hyd is None:
            print("Hydraulic Diameter is not defined")
            return
        if self.k_t is None:
            Re = corr.Re(rho=self.rho, u=self.U0, L=self.d_Hyd, mu=self.mu)
            Sc = corr.Schmidt(D=self.D, mu=self.mu, rho=self.rho)
            if turbulator is None:

                # if Re < 1e4 and Re > 2030:
                #     Sh = 0.015 * Re**0.83 * Sc**0.42  ## Stempien Thesis pg 155-157 TODO implement different Re ranges
                if Re > 2030:
                    Sh = 0.0096 * Re**0.913 * Sc**0.346  ##Getthem paper
                    # Sh = 0.023 * Re**0.8 * Sc**0.33
                else:
                    print(str(Re) + " indicates laminar flow")
                    Sh = 3.66
                    # raise ValueError("Reynolds number is too low")
                self.k_t = corr.get_k_from_Sh(
                    Sh=Sh,
                    L=self.d_Hyd,
                    D=self.D,
                )
            else:
                match turbulator.turbulator_type:
                    case "WireCoil":

                        self.k_t = turbulator.k_t_correlation(
                            Re=Re, Sc=Sc, d_hyd=self.d_Hyd, D=self.D
                        )
                    case "TwistedTape":
                        raise NotImplementedError("Twisted Tape not implemented yet")
                    case "Custom":
                        self.k_t = turbulator.k_t_correlation(
                            Re=Re, Sc=Sc, d_hyd=self.d_Hyd, D=self.D
                        )

        else:
            print("k_t is already defined")


class Membrane(TriomaClass):
    """
    Represents a metallic membrane of a component for H transport.

    Attributes:
        T (float): Temperature of the membrane.
        D (float): Diffusion coefficient of the membrane.
        thick (float): Thickness of the membrane.
        K_S (float): Solubility coefficient of the membrane.
        k_d (float, optional): Dissociation rate constant of the membrane. Defaults to None.
        k_r (float, optional): Recombination rate constant of the membrane. Defaults to None.
        k (float, optional): Thermal conductivity of the membrane. Defaults to None.
        D_0 (float, optional): Pre-exponential factor of the membrane. Defaults to None.Overwrites D if defined
        E_d (float, optional): Activation energy of the diffusivity in the membrane in eV. Defaults to None. Overwrites D if defined
        K_S_0 (float, optional): Pre-exponential factor of the solubility in the membrane. Defaults to None.Overwrites K_S if defined
        E_S (float, optional): Activation energy of the solubility in the membrane in eV. Defaults to None. Overwrites K_S if defined
        inv (float, optional): Inventory of the membrane in mol. Defaults to None.
    """

    def __init__(
        self,
        T: float = None,
        D: float = None,
        thick: float = None,
        K_S: float = None,
        k_d: float = None,
        k_r: float = None,
        k: float = None,
        D_0: float = None,
        E_d: float = None,
        K_S_0: float = None,
        E_S: float = None,
        inv: float = None,
        V: float = None,
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
        if D_0 is not None and E_d is not None:
            self.D = D_0 * np.exp(-E_d / (8.617333262145e-5 * self.T))
        else:
            self.D = D
        if K_S_0 is not None and E_S is not None:
            self.K_S = K_S_0 * np.exp(-E_S / (8.617333262145e-5 * self.T))
        else:
            self.K_S = K_S
        self.thick = thick
        self.k_d = k_d
        self.k_r = k_r
        self.k = k
        self.D_0 = D_0
        self.E_d = E_d
        self.inv = inv
        self.K_S_0 = K_S_0
        self.E_S = E_S

    def set_properties_from_solid_material(
        self, solid_material: "SolidMaterial" = None
    ):
        """
        Sets the properties of the membrane from a SolidMaterial object.

        Args:
            solid_material (SolidMaterial): The SolidMaterial object to set the properties from.
        """
        self.T = solid_material.T
        self.D = solid_material.D
        self.K_S = solid_material.K_S

    def update_T_prop(self):
        if self.D_0 is not None and self.E_d is not None:
            self.D = self.D_0 * np.exp(-self.E_d / (8.617333262145e-5 * self.T))
        if self.K_S_0 is not None and self.E_S is not None:
            self.K_S = self.K_S_0 * np.exp(-self.E_S / (8.617333262145e-5 * self.T))
        return


class FluidMaterial(TriomaClass):
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

    def __init__(
        self,
        T: float = None,
        D: float = None,
        Solubility: float = None,
        MS: bool = None,
        mu: float = None,
        rho: float = None,
        k: float = None,
        cp: float = None,
    ):
        self.T = T
        self.D = D
        self.Solubility = Solubility
        self.MS = MS
        self.mu = mu
        self.rho = rho
        self.k = k
        self.cp = cp


class SolidMaterial(TriomaClass):
    """
    Represents a solid material used in a component.

    Attributes:
        D (float): The Diffusivity of the solid material.
        K_S (float): The Sievert constant of the solid material.
    """

    def __init__(
        self, T: float = None, D: float = None, K_S: float = None, k: float = None
    ):
        self.T = T
        self.D = D
        self.K_S = K_S
        self.k = k


class Turbulator(TriomaClass):
    """
    Represents a turbulator in a component for Tritium transport analysis

    Args:
        turbulator_type (str): Type of the turbulator.
        turbulator_params (dict): Parameters of the turbulator.
    """

    def __init__(
        self,
        turbulator_type: str = None,
    ):
        """
        Initializes a new instance of the Turbulator class.

        Args:
            turbulator_type (str): Type of the turbulator.
            turbulator_params (dict): Parameters of the turbulator.
        """
        self.turbulator_type = turbulator_type


class WireCoil(Turbulator):
    """
    Represents a wire coil turbulator in a component for Tritium transport analysis

    Args:
        turbulator_type (str): Type of the turbulator.
        turbulator_params (dict): Parameters of the turbulator.
    """

    def __init__(
        self,
        turbulator_type: str = "WireCoil",
        pitch: float = None,
    ):
        """
        Initializes a new instance of the WireCoil class.

        Args:
            turbulator_type (str): Type of the turbulator.
            pitch (float): Pitch of the wire coil.
        """
        super().__init__(turbulator_type)
        self.pitch = pitch

    def k_t_correlation(
        self, Re: float = None, Sc: float = None, d_hyd: float = None, D: float = None
    ):
        """
        Calculates the mass transport coefficient (k_t) for the fluid with a wire coil turbulator.

        Args:
            Re (float): Reynolds number of the fluid.
            Sc (float): Schmidt number of the fluid.

        Returns:
            float: The mass transport coefficient.
        """
        if Re > 2030:
            Sh = 0.132 * Re**0.72 * Sc**0.37 * (self.pitch / d_hyd) ** -0.372
        else:
            Sh = 3.66
        k_t = corr.get_k_from_Sh(Sh=Sh, L=self.pitch, D=D)
        return k_t

    def h_t_correlation(
        self, Re: float = None, Pr: float = None, d_hyd: float = None, k=None
    ):
        """
        Calculates the heat transfer coefficient (h_t) for the fluid with a wire coil turbulator.

        Args:
            Re (float): Reynolds number of the fluid.
            Pr (float): Prandtl number of the fluid.

        Returns:
            float: The heat transfer coefficient.
        """
        if Re > 2030:
            Nu = 0.132 * Re**0.72 * Pr**0.37 * (self.pitch / d_hyd) ** -0.372
        else:
            Nu = 3.66
        h_t = corr.get_h_from_Nu(Nu=Nu, k=k, D=d_hyd)
        return h_t


class CustomTurbulator(Turbulator):
    """
    Represents a custom turbulator in a component for Tritium transport analysis
    """

    def __init__(
        self,
        turbulator_type: str = "Custom",
        a: float = None,
        b: float = None,
        c: float = None,
    ):
        """
        Initializes a new instance of the CustomTurbulator class.

        Args:
            turbulator_type (str): Type of the turbulator.
            turbulator_params (dict): Parameters of the turbulator.
        """
        super().__init__(turbulator_type)
        self.a = a
        self.b = b
        self.c = c

    def k_t_correlation(
        self, Re: float = None, Sc: float = None, d_hyd: float = None, D: float = None
    ):
        """
        Calculates the mass transport coefficient (k_t) for the fluid with a custom turbulator.

        Args:
            Re (float): Reynolds number of the fluid.
            Sc (float): Schmidt number of the fluid.

        Returns:
            float: The mass transport coefficient.
        """
        if Re > 2030:
            Sh = self.a * Re**self.b * Sc**self.c
        else:
            Sh = 3.66
        k_t = corr.get_k_from_Sh(Sh=Sh, L=d_hyd, D=D)
        return k_t

    def h_t_correlation(
        self, Re: float = None, Pr: float = None, d_hyd: float = None, k: float = None
    ):
        """
        Calculates the heat transfer coefficient (h_t) for the fluid with a custom turbulator.

        Args:
            Re (float): Reynolds number of the fluid.
            Pr (float): Prandtl number of the fluid.

        Returns:
            float: The heat transfer coefficient.
        """
        if Re > 2030:
            Nu = self.a * Re**self.b * Pr**self.c
        else:
            Nu = 3.66
        h_t = corr.get_h_from_Nu(Nu=Nu, k=k, L=self.pitch)
        return h_t
