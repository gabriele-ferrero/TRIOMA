from TRIOMA.tools.TriomaClass import TriomaClass
from TRIOMA.tools.Extractors.PipeSubclasses import Fluid, Membrane
import TRIOMA.tools.correlations as corr
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
        Calculates the outlet concentration of the GLC.
        The outlet concentration is calculated based on the inlet concentration, the efficiency of the GLC, and the fluid properties.
        Returns:
            None
        """
        match self.fluid.MS:
            case False:
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
                self.eff = eff
                self.c_out = c_out
            case True:
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
                self.eff = eff
                self.c_out = c_out

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
