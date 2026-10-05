import numpy

from TRIOMA.tools.TriomaClass import TriomaClass
from TRIOMA.tools.Extractors.PipeSubclasses import Fluid, Membrane
import TRIOMA.tools.Extractors.extractor as extractor


class GLC_Gas(TriomaClass):
    """
    GLC_Gas class represents a sweep gas in a GLC (Gas-Liquid Contactor) system.
    Attributes:
        G_gas (float): The flow rate of the gas.
        pg_in (float, optional): The Tritium inlet partial pressure of the gas. Default is 0.
        p_tot (float, optional): The total pressure of the component. Default is 100000 Pa.
    """

    def __init__(
        self,
        G_gas: float = None,
        pg_in: float = 0,
        pg_out: float = 0,
        p_tot: float = 100000,
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
        self.pg_out = pg_out


class GLC(TriomaClass):
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
        H: float = None,
        R: float = None,
        L: float = None,
        c_in: float = None,
        eff: float = None,
        fluid: "Fluid" = None,
        membrane: "Membrane" = None,
        GLC_gas: "GLC_Gas" = None,
        T: float = None,
        G_L: float = None,
        c_out: float = None,
        kla: float = None,
    ):
        """
        Initializes a new instance of the GLC class.
        Args:
            H (float): Height of the GLC.
            R (float): Radius of the GLC.
            L (float): Characteristic Length of the GLC fluid flow.
            c_in (float): Inlet concentration of the GLC.
            c_out (float): Outlet concentration of the GLC.
            eff (float, optional): Efficiency of the GLC. Defaults to None.
            fluid (Fluid, optional): Fluid object representing the liquid phase. Defaults to None.
            membrane (Membrane, optional): Membrane object representing the membrane used in the GLC. Not super important. Defaults to None.
            GLC_gas (GLC_Gas, optional): GLC_Gas object representing the gas phase. Defaults to None.
        """
        self.c_in = c_in
        self.eff = eff
        self.fluid = fluid
        self.H = H
        self.R = R
        self.L = L
        self.G_L = G_L
        self.T = T
        self.GLC_gas = GLC_gas
        self.c_out = c_out
        self.eff = eff
        self.kla = kla

    # def get_kla_Ring(self):
    #     """
    #     Calculates the mass transfer coefficient (kla) for a Raschig ring matrix.
    #     The mass transfer coefficient is calculated based on the Reynolds number (Re) and Schmidt number (Sc) of the fluid,
    #     as well as the diameter (d) of the ring.
    #     Returns:
    #         None
    #     """
    #     d = 2e-3  # Ring diameter
    #     Re = corr.Re(rho=self.fluid.rho, u=self.fluid.U0, L=self.L, mu=self.fluid.mu)
    #     Sc = corr.Schmidt(D=self.fluid.D, mu=self.fluid.mu, rho=self.fluid.rho)
    #     self.kla = extractor.corr_packed(
    #         Re,
    #         Sc,
    #         d,
    #         rho_L=self.fluid.rho,
    #         mu_L=self.fluid.mu,
    #         L=self.L,
    #         D=self.fluid.D,
    #     )

    def get_c_out(self):
        """
        Calculate the liquid outlet concentration, liquid extraction
        efficiency, and gas outlet partial pressure.
        """

        if self.fluid is None:
            raise ValueError("A liquid fluid must be assigned to the GLC.")

        if self.GLC_gas is None:
            raise ValueError("A GLC_Gas object must be assigned to the GLC.")

        if self.G_L is None or self.G_L <= 0.0:
            raise ValueError("G_L must be a positive liquid volumetric flow rate.")

        if self.GLC_gas.G_gas is None or self.GLC_gas.G_gas <= 0.0:
            raise ValueError("G_gas must be a positive gas flow rate.")

        if self.H is None or self.H < 0.0:
            raise ValueError("The GLC height H must be non-negative.")

        if self.R is None or self.R <= 0.0:
            raise ValueError("The GLC radius R must be positive.")

        if self.kla is None or self.kla < 0.0:
            raise ValueError("kla must be a non-negative value.")

        R_const = 8.314
        area = numpy.pi * self.R**2

        u_l = self.G_L / area
        u_g = extractor.calculate_gas_velocity(
            G_gas=self.GLC_gas.G_gas,
            p_t=self.GLC_gas.p_tot,
            T=self.T,
            R=self.R,
        )

        if u_g <= 0.0:
            raise ValueError("The calculated gas velocity must be positive.")

        if self.fluid.MS is False:
            # LM concentration follows Sievert's law:
            # c_l = K_S * sqrt(p_l)
            c_out, eff = extractor.get_c_out_GLC_lm(
                Z=self.H,
                R=self.R,
                G_l=self.G_L,
                G_gas=self.GLC_gas.G_gas,
                pl_in=self.c_in**2 / self.fluid.Solubility**2,
                T=self.T,
                p_t=self.GLC_gas.p_tot,
                K_S=self.fluid.Solubility,
                pg_in=self.GLC_gas.pg_in,
                kla=self.kla,
            )

            molecular_factor = 2.0

        else:
            # MS concentration follows Henry's law:
            # c_l = K_H * p_l
            c_out, eff = extractor.get_c_out_GLC_ms(
                Z=self.H,
                R=self.R,
                G_l=self.G_L,
                G_gas=self.GLC_gas.G_gas,
                pl_in=self.c_in / self.fluid.Solubility,
                T=self.T,
                p_t=self.GLC_gas.p_tot,
                K_H=self.fluid.Solubility,
                pg_in=self.GLC_gas.pg_in,
                kla=self.kla,
            )

            molecular_factor = 1.0

        self.eff = float(eff)
        self.c_out = float(c_out)

        # Gas concentration balance:
        #
        # LM: u_g dc_g = -(u_l/2) dc_l
        # MS: u_g dc_g = -u_l dc_l
        #
        # Therefore:
        # c_g,out = c_g,in + u_l/(factor*u_g)*(c_l,in-c_l,out)
        c_g_in = self.GLC_gas.pg_in / (R_const * self.T)

        c_g_out = c_g_in + (u_l / (molecular_factor * u_g)) * (self.c_in - self.c_out)

        self.GLC_gas.pg_out = c_g_out * R_const * self.T

        return self.c_out, self.eff

    def get_kla_from_cout(self):
        match self.fluid.MS:
            case False:
                Bl, kla = extractor.extractor_lm(
                    Z=self.H,
                    R=self.R,
                    G_l=self.G_L,
                    G_gas=self.GLC_gas.G_gas,
                    pl_in=self.c_in**2 / self.fluid.Solubility**2,
                    pl_out=self.c_out**2 / self.fluid.Solubility**2,
                    T=self.T,
                    p_t=self.GLC_gas.p_tot,
                    K_S=self.fluid.Solubility,
                    pg_in=self.GLC_gas.pg_in,
                )
                self.kla = kla
                self.Bl = Bl

            case True:
                Bl, kla = extractor.extractor_ms(
                    Z=self.H,
                    R=self.R,
                    G_l=self.G_L,
                    G_gas=self.GLC_gas.G_gas,
                    pl_in=self.c_in / self.fluid.Solubility,
                    pl_out=self.c_out / self.fluid.Solubility,
                    T=self.T,
                    p_t=self.GLC_gas.p_tot,
                    K_H=self.fluid.Solubility,
                    pg_in=self.GLC_gas.pg_in,
                )
                self.kla = kla
                self.Bl = Bl
        # Update gas outlet partial pressure from the overall isotope balance.
        R_const = 8.314
        area = numpy.pi * self.R**2

        u_l = self.G_L / area
        u_g = extractor.calculate_gas_velocity(
            G_gas=self.GLC_gas.G_gas,
            p_t=self.GLC_gas.p_tot,
            T=self.T,
            R=self.R,
        )

        molecular_factor = 1.0 if self.fluid.MS else 2.0

        c_g_in = self.GLC_gas.pg_in / (R_const * self.T)
        c_g_out = c_g_in + (u_l / (molecular_factor * u_g)) * (self.c_in - self.c_out)

        self.GLC_gas.pg_out = c_g_out * R_const * self.T
        return Bl, kla

    def get_z_from_eff(self):
        """
        Calculates the height of the GLC from the efficiency of the GLC.
        The height is calculated based on the efficiency of the GLC and the radius of the GLC.
        Returns:
            None
        """
        match self.fluid.MS:
            case False:
                z = extractor.length_extractor_lm(
                    R=self.R,
                    G_l=self.G_L,
                    G_gas=self.GLC_gas.G_gas,
                    pl_in=self.c_in**2 / self.fluid.Solubility**2,
                    pl_out=self.c_out**2 / self.fluid.Solubility**2,
                    T=self.T,
                    p_t=self.GLC_gas.p_tot,
                    K_S=self.fluid.Solubility,
                    pg_in=self.GLC_gas.pg_in,
                    kla=self.kla,
                )

            case True:
                z = extractor.length_extractor_ms(
                    R=self.R,
                    G_l=self.G_L,
                    G_gas=self.GLC_gas.G_gas,
                    pl_in=self.c_in / self.fluid.Solubility,
                    pl_out=self.c_out / self.fluid.Solubility,
                    T=self.T,
                    p_t=self.GLC_gas.p_tot,
                    K_H=self.fluid.Solubility,
                    pg_in=self.GLC_gas.pg_in,
                    kla=self.kla,
                )

        return z

    # def connect_to_component(
    #     self, component2: Union["Component", "BreedingBlanket"] = None
    # ):
    #     """sets the inlet conc of the object component equal to the outlet of self"""
    #     component2.update_attribute("c_in", self.c_out)

    def connect_to_component(self):
        return
