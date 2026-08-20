from TRIOMA.tools.Extractors.PipeSubclasses import (
    Geometry,
    Fluid,
    Membrane,
)
from TRIOMA.tools.TriomaClass import TriomaClass
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import lambertw
from scipy import integrate
from typing import Union
from TRIOMA.tools import correlations as corr
import TRIOMA.tools.molten_salts as MS
import TRIOMA.tools.liquid_metals as LM
from scipy.optimize import brentq, least_squares, minimize_scalar


class Component(TriomaClass):
    """
    Represents a component in a plant to make a high level T transport analysis.

    Args:
        Geometry (Geometry): The geometry of the component.
        c_in (float): The concentration of the component at the inlet.
        fluid (Fluid): The fluid associated with the component. Defaults to None.
        membrane (Membrane): The membrane associated with the component. Defaults to None.
    """

    def __init__(
        self,
        geometry: "Geometry" = None,
        c_in: float = None,
        eff: float = None,
        fluid: "Fluid" = None,
        membrane: "Membrane" = None,
        name: str = None,
        p_out: float = 0,
        loss: bool = False,
        inv: float = None,
        delta_p: float = None,
        pumping_power: float = None,
        U: float = None,
        V: float = None,
        cost: float = None,
    ):
        """
        Initializes a new instance of the Component class.

        Args:
            c_in (float): The concentration of the component at the inlet.
            eff (float, optional): The efficiency of the component. Defaults to None.
            L (float, optional): The length of the component. Defaults to None.
            fluid (Fluid, optional): The fluid associated with the component. Defaults to None.
            membrane (Membrane, optional): The membrane associated with the component. Defaults to None.
            name (str, optional): The name of the component. Defaults to None.
            inv (float, optional): The inverse of the efficiency of the component. Defaults to None.
        """
        self.c_in = c_in
        self.geometry = geometry
        self.eff = eff
        self.n_pipes = self.geometry.n_pipes
        self.fluid = fluid
        self.membrane = membrane
        self.name = name
        self.loss = loss
        self.inv = inv
        self.p_out = p_out
        self.delta_p = delta_p
        self.U = U
        self.pumping_power = pumping_power
        self.cost = cost
        self.update_attribute = self.custom_update_attribute

    def custom_update_attribute(self, attr_name: str, new_value: float) -> None:
        """
        Sets the specified attribute to a new value.

        Args:
            attr_name (str): The name of the attribute to set.
            new_value: The new value for the attribute.
        """
        if attr_name == "T":
            if isinstance(self, Component):
                self.fluid.update_attribute(attr_name, new_value)
                self.membrane.update_attribute(attr_name, new_value)
                self.update_T_prop()
                return
            elif hasattr(self, attr_name):
                setattr(self, attr_name, new_value)
                if isinstance(self, Union[Membrane, Fluid]):
                    self.update_T_prop()
                return

        elif hasattr(self, attr_name):
            setattr(self, attr_name, new_value)
            if attr_name == "n_pipes":
                for attr, value in self.__dict__.items():
                    if isinstance(value, object) and hasattr(value, attr_name):

                        setattr(value, attr_name, new_value)
            return
        else:
            for attr, value in self.__dict__.items():
                if isinstance(value, object) and hasattr(value, attr_name):
                    setattr(value, attr_name, new_value)
                    return
        raise ValueError(f"'{attr_name}' is not an attribute of {self.__class__.__name__}")

    def friction_factor(self, Re: float) -> float:
        """
        Calculates the friction factor for the component.

        Args:
            Re (float): Reynolds number.

        Returns:
            float: The friction factor.
        """
        if Re < 2300:
            f = 64 / Re  ## laminar darcy
        else:
            f = 0.316 / Re**0.25  ## Blasius for smooth pipes
        return f

    def update_T_prop(self) -> None:
        """
        Updates the temperature-dependent properties of the fluid and membrane.
        """
        if self.fluid is not None:
            self.fluid.update_T_prop()
        if self.membrane is not None:
            self.membrane.update_T_prop()

    def get_pressure_drop(self) -> float:
        """
        Calculates the pressure drop across the component.

        Returns:
            float: The pressure drop across the component.
        """
        rho = self.fluid.rho
        U = self.fluid.U0
        D = self.geometry.D
        mu = self.fluid.mu
        L = self.geometry.L
        Re = corr.Re(rho, U, D, mu)
        f = self.friction_factor(Re)
        self.delta_p = f * (L / D) * (rho * U**2) / 2
        return self.delta_p

    def estimate_cost(self, metal_cost: float = 0, fluid_cost: float = 0) -> float:
        """
        Estimates the cost of the component.
        metal_cost: cost of the metal in $/m^3
        fluid_cost: cost of the fluid in $/m^3
        returns the cost of the component
        """
        V_solid = self.geometry.get_solid_volume()
        V_fluid = self.geometry.get_fluid_volume()
        cost_solid = V_solid * metal_cost * self.geometry.n_pipes
        cost_fluid = V_fluid * fluid_cost * self.geometry.n_pipes
        self.cost = cost_solid + cost_fluid
        return self.cost

    def get_pumping_power(self) -> float:
        """
        Calculates the pumping power required for the component.

        Returns:
            float: The pumping power required for the component in W.
        """
        if self.delta_p is None:
            self.get_pressure_drop()
        self.pumping_power = self.delta_p * self.get_pipe_flowrate() * self.geometry.n_pipes
        return self.pumping_power

    # def connect_to_component(
    #     self, component2: Union["Component", "BreedingBlanket"] = None
    # ):
    #     """sets the inlet conc of the object component equal to the outlet of self"""
    #     component2.update_attribute("c_in", self.c_out)
    def connect_to_component(self) -> None:
        return  ## empty method defined in component_tools.py

    def plot_component(self) -> plt.Figure:
        r_tot = (self.geometry.D) / 2 + self.geometry.thick
        # Create a figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

        # First subplot: two overlapping circles
        circle1 = plt.Circle((0, 0), r_tot, color="#C0C0C0", label="Membrane")
        circle2 = plt.Circle((0, 0), self.geometry.D / 2, color="#87CEEB", label="Fluid")
        ax1.add_artist(circle1)
        ax1.add_artist(circle2)
        ax1.set_aspect("equal")
        ax1.set_xlim(-r_tot * 1.1, r_tot * 1.1)
        ax1.set_ylim(-r_tot * 1.1, r_tot * 1.1)
        if self.name is None:
            ax1.set_title("Component cross section")
        else:
            ax1.set_title(self.name + " cross section")

        # Add text over the circles
        ax1.text(
            0,
            0,
            r"R=" + str(self.geometry.D / 2),
            color="black",
            ha="center",
            va="center",
        )
        ax1.text(
            0,
            self.geometry.D / 2,
            r"t=" + str(self.geometry.thick),
            color="black",
            ha="center",
            va="center",
            alpha=0.7,
        )
        if self.geometry.n_pipes is not None:
            ax1.text(
                0,
                -self.geometry.D / 4,
                r"" + str(self.geometry.n_pipes) + " pipes",
                color="black",
                ha="center",
                va="center",
                alpha=0.7,
            )

        # Add legend for the circles
        ax1.legend(loc="upper right")
        ax1.axis("off")
        # Second subplot: rectangle and two arrows
        rectangle = plt.Rectangle(
            (0.2, 0.3), 0.55, 0.4, edgecolor="black", facecolor="blue", alpha=0.5
        )
        ax2.add_patch(rectangle)
        # Arrow pointing to the left side of the rectangle
        ax2.arrow(0.0, 0.5, 0.1, 0, head_width=0.05, head_length=0.1, fc="black", ec="black")
        # Arrow pointing out of the right side of the rectangle
        ax2.arrow(0.8, 0.5, 0.1, 0, head_width=0.05, head_length=0.1, fc="black", ec="black")
        ax2.set_aspect("equal")
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        if self.name is None:
            ax2.set_title("Component Lateral view")
        else:
            ax2.set_title(self.name + " Lateral view")

        # Add text over the arrows
        ax2.text(
            0.15,
            0.3,
            r"L=" + str(self.geometry.L) + "m",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.15,
            0.4,
            r"T=" + str(self.fluid.T) + "K",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.15,
            0.6,
            f"c={self.c_in:.4g} $mol/m^3$",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.5,
            0.6,
            f"velocity={self.fluid.U0:.2g} m/s",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.5,
            0.4,
            f"eff={self.eff*100:.2g}%",
            color="black",
            ha="center",
            va="center",
        )

        ax2.text(
            0.9,
            0.3,
            rf"c={self.c_out:.4g}$mol/m^3$",
            color="black",
            ha="center",
            va="center",
        )
        ax2.axis("off")
        # Display the plot
        fig.tight_layout()
        return fig

    def outlet_c_comp(self) -> float:
        """
        Calculate the tritium outlet concentration accounting for extraction and recirculation.

        This method computes the outlet concentration based on the component efficiency and
        inlet concentration, with special handling for feedback effects via recirculation
        (bypass or return flow).

        Returns:
            float: Outlet tritium concentration [mol/m³]
                   Stored in self.c_out

        Three Operating Modes:

            1. **No Recirculation** (recirculation == 0):
                c_out = c_in * (1 - eff)

                Simple extraction with no feedback. The efficiency is applied once.

            2. **Positive Recirculation** (0 < recirculation < 1):
                Solves iteratively for steady-state:
                c_out = c_in * (1 - eff)
                c_in(new) = (c_out * recirculation + c_0) / (recirculation + 1)

                Recirculated tritium in the breeder returns and mixes with fresh inlet stream.
                Uses Picard iteration (tol=1e-6) to reach steady state.

                Example: recirculation=0.5 means 50% of outlet flows back to inlet.
                Physical interpretation: Bypass valve that recycles some extracted tritium

            3. **Bypass without Recirculation** (recirculation < 0, |recirculation| < 1):
                c_out = c_in * (1 - eff) * (1 + recirculation) + c_in * (-recirculation)

                Negative recirculation represents bypass flow: portion of inlet bypasses
                the component entirely and combines with outlet mixture.

                Example: recirculation=-0.3 means 30% of flow bypasses the component
                Physical interpretation: Manifold that splits inlet into two paths

            4. **Invalid Modes**:
                Raises error if recirculation < -1 (bypass exceeds inlet flow)
                Raises error if c_in = 0 and recirculation affects result

        Physics:

            The recirculation parameter models fuel cycle **hydraulic feedback**:

            - **Positive (recycling)**: Reflects scenarios where the processed breeder returns
            in the component and mixes with the fresh inlet stream, increasing inlet concentration
              and potentially improving extraction rate due to higher residence time.

            - **Negative (bypass)**: Represents manifold design where some injected tritium
              bypasses extraction component to improve tritium inventory control

        Parameters Used:
            self.c_in: Inlet concentration [mol/m³]
            self.eff: Component extraction efficiency (dimensionless)
            self.fluid.recirculation: Recirculation coefficient (dimensionless)

        Raises:
            ValueError: If c_in == 0 with active recirculation
            ValueError: If recirculation ≤ -1.0 (bypass exceeds inlet)
            ValueError: If recirculation is NaN or invalid

        Notes:
            - Iteration stops when relative error < 1e-6
            - Maximum iterations: inherent to Picard convergence
            - For recirculation ≠ 0, eff must be pre-calculated (call use_analytical_efficiency first)
        """
        if self.fluid.recirculation == 0:
            self.c_out = self.c_in * (1 - self.eff)
        elif self.fluid.recirculation > 0:
            err = 1
            tol = 1e-6
            if self.c_in == 0:
                raise ValueError("The inlet concentration is zero")
            c0 = self.c_in
            c_in = c0
            while err > tol:

                c_in1 = c_in
                self.update_attribute("c_in", c_in)
                self.analytical_efficiency()
                self.update_attribute("eff", self.eff_an)
                self.c_out = self.c_in * (1 - self.eff)
                c_in = (self.c_out * self.fluid.recirculation + c0) / (self.fluid.recirculation + 1)
                err = abs((c_in - c_in1) / c_in)
        elif self.fluid.recirculation < 0:
            if self.fluid.recirculation <= -1:
                raise ValueError(
                    "Bypass(negative recirculation) not valid: it is more than the flowrate"
                )
            if self.c_in == 0:
                raise ValueError("The inlet concentration is zero")
            self.c_out = self.c_in * (1 - self.eff) * (1 + self.fluid.recirculation) + self.c_in * (
                -self.fluid.recirculation
            )
        else:
            raise ValueError("Recirculation factor not valid")
        return self.c_out

    def converge_split_HX(
        self,
        tol: float = 1e-3,
        T_in_hot: float | None = None,
        T_out_hot: float | None = None,
        T_in_cold: float | None = None,
        T_out_cold: float | None = None,
        R_sec: float | None = None,
        Q: float | None = None,
        plotvar: bool = False,
        savevar: bool = False,
    ) -> None:
        """
        Splits the component into N components to better discretize Temperature effects
        Tries to find the optimal number of components to split the component into

        """

        eff_v = []
        for N in range(10, 101, 2):
            circuit = self.split_HX(
                N=N,
                T_in_hot=T_in_hot,
                T_out_hot=T_out_hot,
                T_in_cold=T_in_cold,
                T_out_cold=T_out_cold,
                R_sec=R_sec,
                Q=Q,
                plotvar=False,
            )

            circuit.get_eff_circuit()
            eff_v.append(circuit.eff)
        x_values = range(10, 101, 2)
        if plotvar == True:
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))

            # First subplot
            axs[0].plot(x_values, eff_v)
            axs[0].set_xlabel("Number of components")
            axs[0].set_ylabel(r"$\eta$")

            # Second subplot
            axs[1].semilogy(
                x_values,
                abs(eff_v - eff_v[-1]) / eff_v[-1] * 100,
                label="Relative error %",
            )
            axs[1].set_xlabel("Number of components")
            axs[1].set_ylabel(r"Relative error (%) in $\eta$ with respect to 100 components")
            axs[1].hlines(
                tol * 100,
                10,
                100,
                colors="r",
                linestyles="dashed",
                label=f"Tolerance {tol*100}%",
            )
            ## remove top and right spines
            axs[0].spines["top"].set_visible(False)
            axs[0].spines["right"].set_visible(False)
            axs[1].spines["top"].set_visible(False)
            axs[1].spines["right"].set_visible(False)
            axs[1].legend(frameon=False, loc="upper right")
            plt.tight_layout()
            if savevar:
                # Save and show the figure
                plt.savefig("HX_convergence.png", dpi=300)
            plt.show()

    def split_HX(
        self,
        N: int = 25,
        T_in_hot: float | None = None,
        T_out_hot: float | None = None,
        T_in_cold: float | None = None,
        T_out_cold: float | None = None,
        R_sec: float = 0,
        Q: float | None = None,
        plotvar: bool = False,
        savevar: bool = False,
    ) -> "Circuit":
        """
        Splits the component into N components to better discretize Temperature effects
        """
        import copy
        from TRIOMA.tools.Circuit import Circuit

        deltaTML = corr.get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold)

        self.get_global_HX_coeff(R_sec)
        ratio_ps = (T_in_hot - T_out_hot) / (
            T_out_cold - T_in_cold
        )  # gets the ratio between flowrate and heat capacity of primary and secondary fluid
        components_list = []

        # Use a for loop to append N instances of Component to the list
        for i in range(N - 1):

            components_list.append(
                Component(
                    name=f"HX_{i+1}",
                    geometry=copy.deepcopy(self.geometry),
                    c_in=copy.deepcopy(self.c_in),
                    fluid=copy.deepcopy(self.fluid),
                    membrane=copy.deepcopy(self.membrane),
                    U=copy.deepcopy(self.U),
                    loss=copy.deepcopy(self.loss),
                )
            )
        T_vec_p = np.linspace(T_in_hot, T_out_hot, N)
        L_vec = []
        position_vec = [0]
        position = 0
        T_vec_s = []
        T_vec_membrane = []
        T_vec_s.append(T_out_cold)
        for i, component in enumerate(components_list):
            deltaTML = corr.get_deltaTML(
                T_in_hot=T_vec_p[i],
                T_out_hot=T_vec_p[i + 1],
                T_in_cold=T_vec_s[i] + (T_vec_p[i + 1] - T_vec_p[i]) / ratio_ps,
                T_out_cold=T_vec_s[i],
            )
            L_vec.append(
                corr.get_length_HX(
                    deltaTML=deltaTML,
                    d_hyd=self.geometry.D,
                    U=component.U,
                    Q=Q / (N - 1),
                )
            )
            next_T_s = T_vec_s[i] + (T_vec_p[i + 1] - T_vec_p[i]) / ratio_ps
            T_vec_s.append(next_T_s)

        for i, component in enumerate(components_list):
            position += L_vec[i]
            position_vec.append(position)
            component.geometry.L = L_vec[i]
            component.fluid.T = (T_vec_p[i] + T_vec_p[i + 1]) / 2
            average_T_s = (T_vec_s[i] + T_vec_s[i + 1]) / 2
            R_prim = 1 / self.fluid.h_coeff
            R_cond = np.log((self.fluid.d_Hyd + self.membrane.thick) / self.fluid.d_Hyd) / (
                2 * np.pi * self.membrane.k
            )
            R_tot = 1 / component.U
            T_membrane = (T_vec_p[i] + T_vec_p[i + 1]) / 2 + (
                average_T_s - ((T_vec_p[i] + T_vec_p[i + 1]) / 2)
            ) * (R_prim + R_cond / 2) / R_tot
            component.membrane.update_attribute("T", T_membrane)
            T_vec_membrane.append(component.membrane.T)
        circuit = Circuit(components=components_list)

        if plotvar:
            fig, axes = plt.subplots(2, 1, figsize=(8, 10))

            # First subplot
            axes[0].plot(T_vec_p)
            axes[0].plot(T_vec_s)
            x_values = np.arange(len(T_vec_membrane)) + 0.5
            axes[0].plot(x_values, T_vec_membrane)
            axes[0].legend(["Primary fluid", "Secondary fluid", "Membrane"], frameon=False)
            axes[0].set_ylabel("Temperature [K]")
            axes[0].set_xlabel("Component number")
            axes[0].spines["top"].set_visible(False)
            axes[0].spines["right"].set_visible(False)

            # Second subplot
            axes[1].plot(position_vec, T_vec_p)
            axes[1].plot(position_vec, T_vec_s)
            axes[1].plot(position_vec[:-1], T_vec_membrane)
            axes[1].legend(["Primary fluid", "Secondary fluid", "Membrane"], frameon=False)
            axes[1].set_xlabel("Position in HX [m]")
            axes[1].set_ylabel("Temperature [K]")
            axes[1].spines["top"].set_visible(False)
            axes[1].spines["right"].set_visible(False)

            # Adjust layout
            plt.tight_layout()
            if savevar:
                # Save and show the figure
                plt.savefig("HX_temperature_profile.png", dpi=300)
            plt.show()

        return circuit

    def T_leak(self) -> float:
        """
        Calculates the leakage of the component.

        Returns:
            float: The leakage of the component.
        """
        leak = self.c_in * self.eff * self.get_pipe_flowrate()

        return leak

    def get_regime(self, print_var: bool = False) -> str:
        """
        Gets the regime of the component.

        Returns:
            str: The regime of the component.
        """
        if self.fluid is None:
            print("No fluid selected")
            return
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        if self.membrane is None:
            print("No membrane selected")
            return
        match self.fluid.MS:
            case True:
                result = MS.get_regime(
                    k_d=self.membrane.k_d,
                    D=self.membrane.D,
                    thick=self.membrane.thick,
                    K_S=self.membrane.K_S,
                    c0=self.c_in,
                    k_t=self.fluid.k_t,
                    k_H=self.fluid.Solubility,
                    print_var=print_var,
                )
                return result
            case False:
                result = LM.get_regime(
                    D=self.membrane.D,
                    k_t=self.fluid.k_t,
                    K_S_S=self.membrane.K_S,
                    K_S_L=self.fluid.Solubility,
                    k_r=self.membrane.k_r,
                    thick=self.membrane.thick,
                    c0=self.c_in,
                    print_var=print_var,
                )
                return result

    def get_pipe_flowrate(self) -> float:
        """
        Calculates the volumetric flow rate of the component [m^3/s].

        Returns:
            float: The flow rate of the component.
        """
        self.pipe_flowrate = self.fluid.U0 * np.pi * self.fluid.d_Hyd**2 / 4
        return self.pipe_flowrate

    def get_total_flowrate(self) -> float:
        """
        Calculates the total flow rate of the component.
        """
        self.get_pipe_flowrate()
        self.flowrate = self.pipe_flowrate * self.geometry.n_pipes
        return self.flowrate

    def define_component_volumes(self) -> None:
        """
        Calculates the volumes of the component.
        """
        self.fluid.V = self.geometry.get_fluid_volume()
        self.membrane.V = self.geometry.get_solid_volume()
        self.V = self.fluid.V + self.membrane.V

    def get_adimensionals(self) -> None:
        """
        Calculate dimensionless transport parameters H and W for tritium permeation analysis.

        These dimensionless numbers characterize the relative importance of different transport
        mechanisms (mass transport vs. diffusion vs. surface kinetics) in tritium permeation.
        They automatically select which physical regime governs extraction and guide flux calculations.

        Dimensionless Parameters:

            **H** (mass transport vs. surface kinetics):
                H = k_t * d_hyd / (k_d * K_S * D)

                - H >> 1: Mass transport is fast → surface reaction becomes rate-limiting
                - H << 1: Surface kinetics are fast → mass transport becomes rate-limiting
                - H ~ 1: Both mechanisms are equally important (mixed regime)

            **W** (diffusion vs. surface kinetics):
                W = (K_S * D / (d_hyd/2)) / k_d  [for molten salts with factor 0.5*K_S*D]

                - W >> 1: Diffusion is slow → surface reaction is fast (diffusion-limited)
                - W << 1: Diffusion is fast → surface reaction is slow (surface-limited)
                - W ~ 1: Both mechanisms coupled (fully mixed regime)

        Fluid Type Corrections:

            **Molten Salt (MS=True)**:
                Uses partition coefficient: K_S = surface/liquid equilibrium
                Includes molecular H₂ dissociation effects in diffusion

            **Liquid Metal (MS=False)**:
                Uses partition coefficient with liquid metal solubility model
                Includes partition parameter for atomic hydrogen transport

        Updates (self attributes):
            self.H (float): Dimensionless parameter (mass transport/surface ratio)
            self.W (float): Dimensionless parameter (diffusion/surface ratio)

        Physics Usage:
            The H and W values automatically route get_flux() calculations to the correct
            transport regime, dramatically reducing computation time:

            If H/W > 1000: Mass transport limited → simple J = -2*k_t*Δc
            If H/W < 0.0001: Diffusion limited → simple J = -(D/δ)*K_S*√(Δc)
            If 0.1 < W < 10: Mixed regime → requires coupled solver

        Dependencies:
            - fluid.k_t: Mass transfer coefficient [m/s] (must be pre-calculated)
            - membrane.k_d: Surface kinetic coefficient [mol/(m²·s)]
            - membrane.K_S: Partition coefficient (dimensionless)
            - membrane.D: Solid-state diffusion coefficient [m²/s]

        Raises:
            None (prints warning if fluid.k_t not yet calculated)
        """
        if self.fluid is None:
            print("No fluid selected")
            return
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        match self.fluid.MS:
            case True:
                self.H = MS.H(k_t=self.fluid.k_t, k_H=self.fluid.Solubility, k_d=self.membrane.k_d)
                self.W = MS.W(
                    k_d=self.membrane.k_d,
                    D=self.membrane.D,
                    thick=self.membrane.thick,
                    K_S=self.membrane.K_S,
                    c0=self.c_in,
                    k_H=self.fluid.Solubility,
                )
            case False:
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

    def use_analytical_efficiency(self, p_out: float = 0) -> None:
        """Evaluates the analytical efficiency and substitutes it in the efficiency attribute of the component.

        Args:
            L (float): the length of the pipe component
        Returns:
            None
        """
        self.analytical_efficiency(p_out=p_out)
        self.eff = self.eff_an

    def get_efficiency(
        self, plotvar: bool = False, c_guess: float | None = None, nodes=100
    ) -> None:
        """
        Calculates the efficiency of the component.
        """
        if self.p_out:
            p_out = self.p_out
        else:
            p_out = 0
        if self.c_in == 0:
            self.c_out = 0
            self.eff = 0
            return

        L_vec = np.linspace(0, self.geometry.L, nodes)
        dl = L_vec[1] - L_vec[0]

        c_vec = np.ndarray(len(L_vec))
        for i in range(len(L_vec)):
            if self.fluid.MS:
                f_H2 = 0.5
            else:
                f_H2 = 1
            if i == 0:

                c_vec[i] = float(self.c_in)

                if isinstance(c_guess, float):
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess, p_out=p_out)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=float(self.c_in), p_out=p_out)
            else:
                c_vec[i] = c_vec[
                    i - 1
                ] + f_H2 * self.J_perm * self.fluid.d_Hyd * np.pi * dl**2 / self.fluid.U0 / (
                    np.pi * self.fluid.d_Hyd**2 / 4 * dl
                )
                if isinstance(c_guess, float):
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess, p_out=p_out)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=float(self.c_in), p_out=p_out)
        if plotvar:
            plt.plot(L_vec, c_vec)
        self.c_out = c_vec[-1]
        self.eff = (self.c_in - self.c_out) / self.c_in

    def analytical_efficiency(self, p_out: float = 0) -> None:
        """
        Calculate the analytical efficiency of a tritium permeation through a component.

        This method computes the tritium extraction efficiency by solving the governing equations
        for tritium transport in the membrane. The efficiency represents the fraction of tritium
        extracted from the component relative to inlet concentration.

        The calculation solves three coupled transport phenomena:
        1. **Mass transport** (fluid boundary layer): Convective mass transfer from bulk fluid to wall
        2. **Diffusion** (solid membrane): Fickian diffusion through the membrane thickness
        3. **Surface reactions** (membrane surfaces): Adsorption/desorption kinetics at interfaces

        Parameters:
            p_out (float): Outlet tritium partial pressure [Pa]. Defaults to 0 Pa (essentially zero).
                           Controls the driving force for tritium extraction.

        Updates (self attributes):
            self.eff_an (float): Analytical efficiency (dimensionless, 0-1)
            self.tau (float): Dimensionless time parameter = 4*k_t*L/(U0*d_Hyd)
            self.alpha (float): Adsorption/surface parameter
            self.xi (float): Extraction parameter

        Physics:
            For **Molten Salt** fluids (MS=True):
                Uses solution of coupled convective-diffusive equations with Lambert W function.
                Handles three limiting regimes: surface-limited, diffusion-limited, and mass-transport-limited.

            For **Liquid Metal** fluids (MS=False):
                Uses simplified solution based on partition equilibrium effects.
                Includes pressure correction factor: (1 - p_out/p_in)^0.5

        References:
            Humrickhouse, P. W., "Tritium Transport in the DCLL Blanket",
            18th ANS Topical Meeting on Fusion Energy, 2008.

        Raises:
            ValueError: If imaginary component appears in eff_an calculation (numerical instability)
        """
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        self.tau = 4 * self.fluid.k_t * self.geometry.L / (self.fluid.U0 * self.fluid.d_Hyd)
        match self.fluid.MS:
            case True:  # Molten salt

                KH = self.fluid.Solubility
                d = self.fluid.d_Hyd
                r_i = d / 2.0
                r_o = r_i + self.membrane.thick
                log_ro_ri = np.log(r_o / r_i)

                phi = self.membrane.D * self.membrane.K_S

                self.alpha = 1.0 / KH * (phi / (self.fluid.k_t * d * log_ro_ri)) ** 2

                self.xi = self.alpha / self.c_in

                if p_out < 0:
                    raise ValueError("p_out must be non-negative.")

                p_in = self.c_in / KH

                self.Pi_ext = np.sqrt(p_out * KH / self.alpha)

                # -------------------------------------------------------------
                # Mass-transfer-limited approximation: xi >> 1
                # -------------------------------------------------------------
                if self.xi > 1.0e5:
                    correction_p = 1.0 - p_out / p_in

                    self.eff_an = (1.0 - np.exp(-self.tau)) * correction_p

                # -------------------------------------------------------------
                # Diffusion-limited approximation:
                # xi << 1 and tau < 1 / sqrt(xi)
                # -------------------------------------------------------------
                elif self.xi < 1.0e-4 and self.tau < 1.0 / np.sqrt(self.xi):
                    correction_p = 1.0 - np.sqrt(p_out / p_in)

                    self.eff_an = (
                        1.0 - (1.0 - 0.5 * self.tau * np.sqrt(self.xi)) ** 2
                    ) * correction_p

                # -------------------------------------------------------------
                # Exact Lambert-W solution
                # -------------------------------------------------------------
                else:
                    Pi_ext = self.Pi_ext
                    b = 1.0 + 2.0 * Pi_ext

                    s_in = np.sqrt(1.0 + 4.0 * (1.0 / self.xi + Pi_ext))

                    y_in = s_in - b

                    # No driving force: p_out = p_in
                    if abs(y_in) < 1.0e-14:
                        self.eff_an = 0.0
                        return

                    beta = s_in / b + np.log(abs(y_in))

                    beta_tau = beta - self.tau / b - 1.0
                    # Lambert-W argument:
                    # extraction:     exp(beta_tau) / b
                    # inverse permeation:
                    #                 -exp(beta_tau) / b
                    log_argument_abs = beta_tau - np.log(b)
                    if y_in > 0.0:
                        # Normal extraction: principal branch W_0
                        log_max = np.log(np.finfo(np.float64).max)
                        if log_argument_abs < log_max:
                            argument = np.exp(log_argument_abs)
                            q_out = lambertw(argument, k=0).real
                        else:
                            # Large-positive-argument asymptotic
                            q_out = log_argument_abs - np.log(log_argument_abs)
                    else:
                        # Inverse permeation:
                        # the physical solution remains on W_0
                        argument = -np.exp(log_argument_abs)
                        # Protect against round-off below -1/e
                        argument = np.clip(argument, -1.0 / np.e, 0.0)
                        q_out = lambertw(argument, k=0).real
                    s_out = b * (1.0 + q_out)
                    c_out_over_alpha = (s_out**2 - 1.0 - 4.0 * Pi_ext) / 4.0

                    self.eff_an = 1.0 - self.xi * c_out_over_alpha

                    return
                    # e = (self.alpha * p_out * self.fluid.Solubility) ** 0.5
                    # f = e / self.alpha
                    # delta = (1 / self.xi + 1 + 2 * f) ** 0.5
                    # beta = delta + (1 + f) * np.log(abs(delta - 1 - f))
                    # print("beta is ", beta)
                    # max_exp = np.log(np.finfo(np.float64).max)
                    # beta_tau = beta - self.tau - 1
                    # print("beta tau is ", beta_tau)
                    # print("max exp is ", max_exp)
                    # saturation = (p_out * self.fluid.Solubility) > self.c_in
                    # if beta_tau > max_exp :
                    #     # we can use the approximation w=beta_tau-np.log(beta_tau)for the lambert W function but it leads to error up to 40 % in very niche scenarios.

                    #     def eq(var):
                    #         cl = var
                    #         self.alpha = self.xi * self.c_in

                    #         left = (cl / self.alpha + 1 + 2 * f) ** 0.5 + (1 + f) * np.log(
                    #             abs(-f + ((cl/self.alpha + 1 + 2*f)**0.5 - 1))
                    #         )

                    #         right = beta - self.tau

                    #         return abs(left - right)

                    #     p_in = self.c_in / self.fluid.Solubility
                    #     if (
                    #         abs(self.p_out * self.fluid.Solubility - self.c_in) / self.c_in
                    #         < 1e-2
                    #     ):
                    #         self.eff_an = 1e-6
                    #         return
                    #     lower_bound = min(self.p_out * self.fluid.Solubility, self.c_in)
                    #     upper_bound = max(self.p_out * self.fluid.Solubility, self.c_in)
                    #     cl = minimize(
                    #         eq,
                    #         x0=(lower_bound + upper_bound) / 2,
                    #         method="Powell",
                    #         bounds=[(lower_bound, upper_bound)],
                    #         tol=1e-7,
                    #     ).x[0]
                    #     # corr_p=1-(p_out/p_in)
                    #     self.eff_an = 1 - (cl / self.c_in)
                    #     return
                    # else:
                    #     z = np.exp(beta_tau)
                    #     w = lambertw(z, tol=1e-10)
                    #     self.eff_an = 1 - self.xi * (w**2 + 2 * w)
                    #     if self.eff_an.imag != 0:
                    #         raise ValueError("self.eff_an has a non-zero imaginary part")
                    #     else:
                    #         self.eff_an = self.eff_an.real  # get rid of 0*j
            case False:  # Liquid Metal
                self.zeta = (2 * self.membrane.K_S * self.membrane.D) / (
                    self.fluid.k_t
                    * self.fluid.Solubility
                    * self.fluid.d_Hyd
                    * np.log((self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd)
                )
                p_in = (self.c_in / self.fluid.Solubility) ** 2
                corr_p = 1 - (p_out / p_in) ** 0.5

                self.eff_an = (1 - np.exp(-self.tau * self.zeta / (1 + self.zeta))) * corr_p

    def _diff_conductance(self) -> float:
        """Cylindrical diffusion conductance D*K_S / (r ln((r+t)/r)), computed once."""
        r = self.fluid.d_Hyd / 2.0
        return self.membrane.D * self.membrane.K_S / (r * np.log((r + self.membrane.thick) / r))

    def get_flux(self, c: float | None = None, c_guess: float = 1e-9, p_out: float = 0) -> float:
        """
        Tritium permeation flux across the membrane.

        Solves for the wall/interface concentration satisfying the coupled
        mass-transport / diffusion / surface-reaction fluxes, stores the result
        in ``self.J_perm`` [mol/(m^2 s)] and returns the wall concentration.



        Raises
        ------
        ValueError
            If ``c`` or ``c_guess`` is not a float.
        """
        if not isinstance(c, float):
            print(c)
            raise ValueError("Input 'c' must be a float")
        if not isinstance(c_guess, float):
            raise ValueError("c_guess must be a float")

        self.get_adimensionals()

        # --- constants -------------------------------------------------------
        S = self.fluid.Solubility
        kt = self.fluid.k_t
        kd = self.membrane.k_d
        KS = self.membrane.K_S
        P = self._diff_conductance()  # D*K_S / (r ln((r+t)/r))

        if self.fluid.MS:  # H2 <-> 2H, Henry
            mt = 2.0
            exp = 0.5  # (cw/S)**0.5 in diffusion driving force
            des_pow = 2.0  # desorption uses K_S**2 (fully-mixed)
            c_eq = S * p_out  # wall conc. in equilibrium with p_out
            lo, hi = min(c, c_eq), max(c, c_eq)
        else:  # atomic H, Sieverts
            mt = 1.0
            exp = 1.0  # (cw/S) in diffusion driving force
            des_pow = 1.0  # desorption uses K_S (fully-mixed)
            c_eq = S * np.sqrt(p_out)
            lo = min(c, c_eq)
            hi = max(c, c_eq)

        # --- signed fluxes (positive = leaving the fluid) --------------------
        def J_mt(cw):
            return mt * kt * (c - cw)

        def J_diff(cw):
            return P * ((cw / S) ** exp - p_out**0.5)

        def J_surf(cw):
            return kd * (c / S) - kd * KS**2 * cw**2

        def root(f, g):
            """cw where f==g on [lo, hi]. brentq if bracketed, else bounded |f-g| min."""
            if hi - lo <= 1e-12 * max(abs(hi), 1e-30):  # no driving force
                return c
            h = lambda cw: f(cw) - g(cw)
            flo, fhi = h(lo), h(hi)
            if flo == 0.0:
                return lo
            if fhi == 0.0:
                return hi
            if flo * fhi < 0.0:  # sign change -> exact root
                return brentq(h, lo, hi, xtol=1e-18, rtol=1e-12, maxiter=200)
            res = minimize_scalar(  # no bracket -> minimise |residual|
                lambda cw: abs(h(cw)),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-14},
            )
            return float(res.x)

        # =====================================================================
        #  W > 10 : diffusion vs mass transport (surface fast)
        # =====================================================================
        if self.W > 10:
            if self.H / self.W > 1000:  # mass-transport limited
                cw = c_eq
                self.J_perm = -J_mt(cw)  # NEGATIVE (leaving the fluid)
            elif self.H / self.W < 1e-4:  # diffusion limited
                cw = c
                self.J_perm = -J_diff(cw)  # NEGATIVE (leaving the fluid)
            else:  # mixed MT <-> diffusion
                cw = root(J_mt, J_diff)
                self.J_perm = -J_mt(cw)  # NEGATIVE (leaving the fluid)
            return float(cw)

        # =====================================================================
        #  W < 0.1 : surface reaction vs mass transport (diffusion fast)
        # =====================================================================
        if self.W < 0.1:
            if self.H > 100:  # mass-transport limited
                cw = c_eq
                self.J_perm = -J_mt(cw)  # NEGATIVE (leaving the fluid)
            elif self.H < 1e-2:  # surface limited
                cw = c
                self.J_perm = -kd * (c / S)  # NEGATIVE (leaving the fluid)
            else:  # mixed MT <-> surface
                cw = root(J_mt, J_surf)
                self.J_perm = -J_mt(cw)  # negative (leaving the fluid)
            return float(cw)

        # =====================================================================
        #  Intermediate W : coupled surface + diffusion (mass transport fast)
        # =====================================================================
        if self.H / self.W > 1000:  # mass-transport limited
            cw = c_eq
            self.J_perm = -J_mt(cw)  # NEGATIVE (leaving the fluid)
            return float(cw)

        if self.H / self.W < 1e-4:  # surface <-> diffusion
            cw = root(J_surf, J_diff)
            self.J_perm = -J_diff(cw)  # NEGATIVE (leaving the fluid)
            return float(cw)

        # --- fully coupled: mass transport + surface + diffusion -------------
        # Desorption K_S power differs by fluid (MS: K_S**2, LM: K_S) -> des_pow.
        def system(v):
            cw, cs = v
            Jmt = mt * kt * (c - cw)
            Jd = kd * (cw / KS) - kd * (KS**des_pow) * cs**2
            Jdiff = P / KS * (cs - KS * p_out**0.5)  # P already contains K_S
            return [Jmt - Jd, Jmt - Jdiff]

        if hi - lo <= 1e-12 * max(abs(hi), 1e-30):
            self.J_perm = 0.0
            return float(c)

        cw0 = float(np.clip(2.0 * c / 3.0, lo, hi))  # original guess, clamped
        cs0 = max(c / 3.0, 0.0)
        sol = least_squares(
            system,
            [cw0, cs0],
            bounds=([lo, 0.0], [hi, np.inf]),
            xtol=1e-12,
            ftol=1e-12,
            max_nfev=int(1e4),
        )
        cw = sol.x[0]
        self.J_perm = -J_mt(cw)  # NEGATIVE (leaving the fluid)
        return float(cw)

    def get_global_HX_coeff(self, R_conv_sec: float = 0) -> None:
        """
        Calculate the overall heat transfer coefficient for a heat exchanger component.

        This method computes the global heat transfer coefficient (U-value) accounting for
        all thermal resistances in series: primary-side convection, membrane conduction,
        and optional secondary-side convection.

        The overall heat transfer is modeled as thermal resistors in series:
        U = 1 / (R_conv_prim + R_cond + R_conv_sec)

        Parameters:
            R_conv_sec (float): Secondary-side convection thermal resistance [K/W].
                               Default 0 (adiabatic or negligible resistance).
                               Represents heat transfer resistance on downstream side.

        Calculates:
            1. **Primary-side convection resistance** R_conv_prim:
               - Determines Nusselt number via appropriate correlation:
                 * Dittus-Boelert (smooth pipes): Nu = 0.023*Re^0.8*Pr^0.4
                 * WireCoil turbulator: custom correlation (if installed)
                 * CustomTurbulator: user-defined correlation
               - Converts Nu to convection coefficient: h = Nu*k/d_hyd
               - R_conv_prim = 1/h

            2. **Membrane conduction resistance** R_cond:
               - Cylindrical geometry: R_cond = ln(r_outer/r_inner) / (2π*k)
               - k: membrane thermal conductivity [W/(m·K)]
               - r_outer/r_inner: outer/inner radii including thickness

        Parameters Used:
            self.fluid: FluidMaterial with properties (ρ, μ, k, cp for correlations)
            self.geometry: Component geometry (D, L, turbulator type)
            self.membrane: SolidMaterial with thermal conductivity k

        Updates (self attributes):
            self.U (float): Overall HX coefficient [W/(m²·K)]
            self.fluid.h_coeff (float): Primary convection coefficient [W/(m²·K)]

        Physics Correlations:
            **Reynolds number**: Re = ρ*U*d_hyd/μ  (flow regime indicator)
            **Prandtl number**: Pr = cp*μ/k  (thermal property ratio)
            **Nusselt number**: dimensionless heat transfer (depends on Re, Pr, geometry)

        Physics/Engineering Note:
            This U-value is used in heat exchanger finite-difference splitting (split_HX)
            to discretize temperature profiles and improve tritium extraction efficiency
            calculations that depend on local temperatures.

        Raises:
            NotImplementedError: If turbulator_type is "TwistedTape" (not yet implemented)
        """
        R_cond = np.log((self.fluid.d_Hyd + self.membrane.thick) / self.fluid.d_Hyd) / (
            2 * np.pi * self.membrane.k
        )
        Re = corr.Re(self.fluid.rho, self.fluid.U0, self.fluid.d_Hyd, self.fluid.mu)
        Pr = corr.Pr(self.fluid.cp, self.fluid.mu, self.fluid.k)
        if self.geometry.turbulator is None:
            h_prim = corr.get_h_from_Nu(
                corr.Nu_DittusBoelter(Re, Pr), self.fluid.k, self.fluid.d_Hyd
            )
        else:
            match self.geometry.turbulator.turbulator_type:
                case "TwistedTape":
                    print(str(self.geometry.turbulator.turbulator_type) + " is not implemented yet")
                    raise NotImplementedError("Twisted tape is not implemented yet")
                case "WireCoil":
                    h_prim = self.geometry.turbulator.h_t_correlation(
                        Re=Re, Pr=Pr, d_hyd=self.fluid.d_Hyd, k=self.fluid.k
                    )
                case "Custom":
                    h_prim = self.geometry.turbulator.h_t_correlation(
                        Re=Re, Pr=Pr, d_hyd=self.fluid.d_Hyd, k=self.fluid.k
                    )
        self.fluid.h_coeff = h_prim
        R_conv_prim = 1 / h_prim
        R_tot = R_conv_prim + R_cond + R_conv_sec
        self.U = 1 / R_tot
        return

    def analytical_solid_inventory(self, p_out: float = 0) -> float:
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        match self.fluid.MS:
            case False:

                def integralfun(r):
                    return (
                        1
                        / 4
                        * r**2
                        * (2 * np.log(r / (self.geometry.D / 2 + self.geometry.thick)) - 1)
                    )

                def circle(r):
                    return np.pi * r**2

                dimless = (
                    2
                    * self.membrane.D
                    * self.membrane.K_S
                    / (
                        self.fluid.k_t
                        * self.fluid.Solubility
                        * self.fluid.d_Hyd
                        * np.log((self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd)
                    )
                )
                dimless2 = (
                    2
                    * self.membrane.D
                    * self.membrane.K_S
                    / (
                        self.fluid.Solubility
                        * self.fluid.d_Hyd
                        * np.log((self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd)
                    )
                )
                L_ch = (
                    -dimless
                    / (1 + dimless)
                    * 4
                    * self.fluid.k_t
                    / (self.fluid.U0 * self.fluid.d_Hyd)
                )
                L_factor = (np.exp(L_ch * self.geometry.L) - 1) / L_ch
                K = (
                    -2
                    * np.pi
                    * (
                        (self.c_in - self.fluid.Solubility * p_out**0.5)
                        / (dimless2 / self.fluid.k_t + 1)
                        / self.fluid.Solubility
                        * self.membrane.K_S
                    )
                    / np.log((self.geometry.D / 2 + self.geometry.thick) / (self.geometry.D / 2))
                ) * L_factor
                integral = (
                    K * integralfun(self.geometry.D / 2 + self.geometry.thick)
                    - K * integralfun(self.geometry.D / 2)
                ) + self.geometry.L * p_out**0.5 * self.membrane.K_S * (
                    circle(self.geometry.D / 2 + self.geometry.thick) - circle(self.geometry.D / 2)
                )
                inventory = integral
                self.membrane.inv = inventory * self.geometry.n_pipes
                return inventory
            case True:

                def ms_integral(self, p_out: float = 0.0, L: float = None):
                    """
                    Solid MS inventory for one pipe, using the paper's
                    alpha, xi and Pi_ext definitions.
                    """

                    if L is None:
                        L = self.geometry.L

                    if not (
                        hasattr(self, "alpha")
                        and hasattr(self, "xi")
                        and self.alpha is not None
                        and self.xi is not None
                    ):
                        self.analytical_efficiency(p_out=p_out)

                    KH = self.fluid.Solubility
                    KS = self.membrane.K_S
                    kt = self.fluid.k_t
                    U = self.fluid.U0
                    d = self.fluid.d_Hyd

                    r_i = d / 2.0
                    r_o = r_i + self.membrane.thick
                    log_ro_ri = np.log(r_o / r_i)

                    alpha = self.alpha
                    xi = self.xi

                    # Paper definition:
                    # Pi_ext = sqrt(p_out * K_H / alpha)
                    Pi_ext = np.sqrt(p_out * KH / alpha)
                    b = 1.0 + 2.0 * Pi_ext

                    # Initial transformed variable
                    s_in = np.sqrt(1.0 + 4.0 * (1.0 / xi + Pi_ext))
                    y_in = s_in - b

                    if abs(y_in) < 1.0e-14:
                        # No concentration driving force
                        c_w_s_integral = KS * np.sqrt(p_out) * L
                    else:
                        sign = 1.0 if y_in > 0.0 else -1.0

                        # Paper definition of beta
                        beta = s_in / b + np.log(abs(y_in))

                        def q_of_z(z):
                            tau_z = 4.0 * kt * z / (U * d)

                            beta_z = beta - tau_z / b - 1.0

                            # q = W_k[sign * exp(beta_z) / b]
                            log_argument_abs = beta_z - np.log(b)

                            if sign > 0.0:
                                # Normal extraction: W_0
                                log_max = np.log(np.finfo(np.float64).max)

                                if log_argument_abs < log_max:
                                    argument = np.exp(log_argument_abs)
                                    q = lambertw(argument, k=0).real
                                else:
                                    # Large-positive-argument
                                    # asymptotic approximation
                                    q = log_argument_abs - np.log(log_argument_abs)

                            else:
                                # Inverse permeation: physical branch W_0
                                if log_argument_abs < np.log(np.finfo(np.float64).tiny):
                                    argument = 0.0
                                else:
                                    argument = -np.exp(log_argument_abs)

                                argument = np.clip(argument, -1.0 / np.e, 0.0)

                                q = lambertw(argument, k=0).real

                            return q

                        def c_w_l(z):
                            """
                            MS liquid-side wall concentration:

                            c_w,l = alpha *
                                    [Pi_ext +
                                     (1 + 2 Pi_ext) q / 2]^2
                            """
                            q = q_of_z(z)

                            return alpha * (Pi_ext + 0.5 * b * q) ** 2

                        def c_w_s(z):
                            """
                            MS solid-side wall concentration:

                            c_w,s = K_S * sqrt(c_w,l / K_H)
                            """
                            cwl = max(c_w_l(z), 0.0)

                            return KS * np.sqrt(cwl / KH)

                        # Numerical axial integration of c_w,s(z)
                        c_w_s_integral, _ = integrate.quad(
                            c_w_s,
                            0.0,
                            L,
                            epsabs=1.0e-12,
                            epsrel=1.0e-8,
                            limit=200,
                        )

                    # Radial hollow-cylinder geometry
                    area_solid = np.pi * (r_o**2 - r_i**2)

                    F_cyl = np.pi * ((r_o**2 - r_i**2) / (2.0 * log_ro_ri) - r_i**2)

                    # External-equilibrium concentration in the solid
                    c_ext_s = KS * np.sqrt(p_out)

                    # Solid inventory in one pipe:
                    #
                    # I_s = c_ext,s * A_s * L
                    #       + F_cyl * integral[c_w,s(z)-c_ext,s] dz
                    inventory_one_pipe = c_ext_s * area_solid * L + F_cyl * (
                        c_w_s_integral - c_ext_s * L
                    )
                    self.membrane.inv = inventory_one_pipe * self.geometry.n_pipes
                    return self.membrane.inv

                inv = ms_integral(
                    self=self,
                    p_out=p_out,
                    L=self.geometry.L,
                )

                self.membrane.inv = inv * self.geometry.n_pipes

                return inv

    def get_solid_inventory(self, p_out: float = 0, flag_an: bool = False) -> float:
        if flag_an:
            return self.analytical_solid_inventory(p_out=p_out)

        def integrate_c_profile(self):
            r_in = self.fluid.d_Hyd / 2
            r_out = self.fluid.d_Hyd / 2 + self.membrane.thick
            L_min = 0
            L_max = self.geometry.L
            N = 20

            def integrand(r, L, p_out=p_out):
                # return -c * np.log(r / r_out) / np.log(r_out / r_in) * 2 * np.pi * r
                if self.fluid.k_t is None:

                    self.fluid.get_kt(turbulator=self.geometry.turbulator)
                if self.fluid.MS == False:
                    c = self.c_in / self.fluid.Solubility * self.membrane.K_S
                    dimless = (
                        2
                        * self.membrane.D
                        * self.membrane.K_S
                        / (
                            self.fluid.k_t
                            * self.fluid.Solubility
                            * self.fluid.d_Hyd
                            * np.log(
                                (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                            )
                        )
                    )
                    dimless2 = (
                        2
                        * self.membrane.D
                        * self.membrane.K_S
                        / (
                            self.fluid.Solubility
                            * self.fluid.d_Hyd
                            * np.log(
                                (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                            )
                        )
                    )
                    L_ch = (
                        -dimless
                        / (1 + dimless)
                        * 4
                        * self.fluid.k_t
                        / (self.fluid.U0 * self.fluid.d_Hyd)
                    )
                    conv_liquid_to_solid = self.membrane.K_S / self.fluid.Solubility
                    c_ext = p_out**0.5 * self.membrane.K_S
                    c_w = (
                        c * np.exp(L_ch * L) / (dimless2 / self.fluid.k_t + 1) + c_ext
                    )  # todo check this is liquid conc

                    return (
                        (-(c_w - c_ext) * np.log(r / r_out) / np.log(r_out / r_in) + c_ext)
                        * 2
                        * np.pi
                        * r
                    )
                else:
                    tau = 4 * self.fluid.k_t * L / (self.fluid.U0 * self.fluid.d_Hyd)
                    self.xi = (
                        1
                        / self.c_in
                        / self.fluid.Solubility
                        * (
                            0.5  ##TODO: Check this
                            * self.membrane.K_S
                            * self.membrane.D
                            / (
                                self.fluid.k_t
                                * self.fluid.d_Hyd
                                * np.log(
                                    (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                                )
                            )
                        )
                        ** 2
                    )

                    beta = (1 / self.xi + 1) ** 0.5 + np.log((1 / self.xi + 1) ** 0.5 - 1)
                    max_exp = np.log(np.finfo(np.float64).max)
                    beta_tau = beta - tau - 1
                    if beta_tau > max_exp:

                        w = beta_tau - np.log(beta_tau)

                    else:
                        z = np.exp(beta_tau)
                        w = lambertw(z, tol=1e-10)
                        if w.imag != 0:
                            raise ValueError("self.eff_an has a non-zero imaginary part")
                        w = w.real
                    alpha = (
                        1
                        / self.fluid.Solubility
                        * (
                            (0.5 * self.membrane.D * self.membrane.K_S)  ## TODO: Check this
                            / (
                                self.fluid.k_t
                                * self.fluid.d_Hyd
                                * np.log(
                                    (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                                )
                            )
                        )
                        ** 2
                    )
                    c_ext = p_out**0.5 * self.membrane.K_S
                    conv = (self.c_in / self.fluid.Solubility) ** 0.5 * self.membrane.K_S
                    c_w_l = (
                        alpha * (w**2 + 2 * w)
                        + alpha * (2 - 2 * ((w**2 + 2 * w) + 1) ** 0.5)  ## TODO: Check this
                        + c_ext
                    )

                    if c_w_l < 0:
                        c_w_l = 1e-17
                    return (
                        (
                            -np.log(r / r_out)
                            / np.log(r_out / r_in)
                            * (
                                (alpha / self.fluid.Solubility) ** 0.5 * w * self.membrane.K_S
                                - c_ext
                            )
                            + c_ext
                        )
                        * 2
                        * np.pi
                        * r
                    )

            result, err = integrate.nquad(integrand, [[r_in, r_out], [L_min, L_max]])
            return result

        integral_pipe = integrate_c_profile(self)
        self.membrane.inv = integral_pipe * self.geometry.n_pipes
        if math.isnan(self.membrane.inv):
            print("Error: Inventory calculation failed")
            self.inspect()
        return self.membrane.inv

    def analytical_fluid_inventory(self, p_out: float = 0) -> None:
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        match self.fluid.MS:
            case False:

                def circle(r):
                    return np.pi * r**2

                dimless = (
                    2
                    * self.membrane.D
                    * self.membrane.K_S
                    / (
                        self.fluid.k_t
                        * self.fluid.Solubility
                        * self.fluid.d_Hyd
                        * np.log((self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd)
                    )
                )
                dimless2 = (
                    2
                    * self.membrane.D
                    * self.membrane.K_S
                    / (
                        self.fluid.Solubility
                        * self.fluid.d_Hyd
                        * np.log((self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd)
                    )
                )
                L_ch = (
                    -dimless
                    / (1 + dimless)
                    * 4
                    * self.fluid.k_t
                    / (self.fluid.U0 * self.fluid.d_Hyd)
                )
                c_ext = p_out**0.5 * self.fluid.Solubility
                L_factor = (np.exp(L_ch * self.geometry.L) - 1) / L_ch
                integral = (self.c_in - c_ext) * circle(
                    self.geometry.D / 2
                ) * L_factor + c_ext * circle(self.geometry.D / 2) * self.geometry.L
                inventory = integral
                self.fluid.inv = inventory * self.geometry.n_pipes
                return inventory
            case True:
                print("MS fluid integration is done numerically")
                self.get_fluid_inventory(flag_an=False, p_out=p_out)

    def get_fluid_inventory(self, flag_an: bool = False, p_out: float = 0) -> float:
        """
        Calculate the tritium inventory in the fluid region.

        For liquid metals, the inventory is calculated analytically from the
        exponential axial concentration profile.

        For molten salts, the bulk concentration is evaluated using the same
        Lambert-W solution used by analytical_efficiency(), then integrated
        numerically along the pipe.

        The molten-salt inventory is multiplied by 2 to convert from mol Q2
        to mol Q, consistently with the manuscript's f_H_to_H2 factor.
        """

        if flag_an:
            return self.analytical_fluid_inventory(p_out=p_out)

        if p_out < 0:
            raise ValueError("p_out must be non-negative.")

        if self.fluid.k_t is None:
            self.fluid.get_kt(turbulator=self.geometry.turbulator)

        c_in = float(self.c_in)
        kt = float(self.fluid.k_t)
        U = float(self.fluid.U0)
        d = float(self.fluid.d_Hyd)
        r_in = d / 2.0
        area_fluid = np.pi * r_in**2
        L = float(self.geometry.L)
        n_pipes = self.geometry.n_pipes

        if c_in < 0:
            raise ValueError("The inlet concentration must be non-negative.")

        # ---------------------------------------------------------------------
        # Liquid-metal carrier: Sievert's law
        # ---------------------------------------------------------------------
        if not self.fluid.MS:

            K_S_l = float(self.fluid.Solubility)
            K_S_s = float(self.membrane.K_S)
            D_s = float(self.membrane.D)
            thickness = float(self.membrane.thick)

            log_ro_ri = np.log((d + 2.0 * thickness) / d)

            # Manuscript definition of zeta
            zeta = 2.0 * D_s * K_S_s / (d * log_ro_ri * kt * K_S_l)

            # Axial decay coefficient:
            #
            # c_b(z) = c_ext + (c_in - c_ext) exp(-a*z)
            #
            a = 4.0 * kt / (U * d) * zeta / (1.0 + zeta)

            c_ext = K_S_l * np.sqrt(p_out)

            if abs(a) < 1.0e-14:
                axial_integral = c_in * L
            else:
                axial_integral = c_ext * L + (c_in - c_ext) * (1.0 - np.exp(-a * L)) / a

            inventory_one_pipe = area_fluid * axial_integral

            self.fluid.inv = inventory_one_pipe * n_pipes
            return self.fluid.inv

        # ---------------------------------------------------------------------
        # Molten-salt carrier: Henry's law and Lambert-W solution
        # ---------------------------------------------------------------------

        if c_in == 0.0:
            if p_out == 0.0:
                self.fluid.inv = 0.0
                return 0.0
            raise ValueError(
                "A positive c_in is required for the molten-salt Lambert-W "
                "inventory formulation when p_out > 0."
            )

        K_H = float(self.fluid.Solubility)
        K_S = float(self.membrane.K_S)
        D_s = float(self.membrane.D)

        r_o = r_in + float(self.membrane.thick)
        log_ro_ri = np.log(r_o / r_in)

        # These definitions must match analytical_efficiency()
        phi = D_s * K_S

        alpha = 1.0 / K_H * (phi / (kt * d * log_ro_ri)) ** 2

        xi = alpha / c_in
        Pi_ext = np.sqrt(p_out * K_H / alpha)

        # Store the parameters for consistency with the rest of the class
        self.alpha = alpha
        self.xi = xi
        self.Pi_ext = Pi_ext
        self.tau = 4.0 * kt * L / (U * d)

        b = 1.0 + 2.0 * Pi_ext

        # Initial transformed variable
        s_in = np.sqrt(1.0 + 4.0 * (1.0 / xi + Pi_ext))
        y_in = s_in - b

        # If the inlet is already in equilibrium with the external pressure,
        # the bulk concentration is constant along the pipe.
        if abs(y_in) < 1.0e-14:
            axial_integral = c_in * L

        else:
            sign = 1.0 if y_in > 0.0 else -1.0

            # Same beta definition as analytical_efficiency()
            beta = s_in / b + np.log(abs(y_in))

            log_max = np.log(np.finfo(np.float64).max)
            log_tiny = np.log(np.finfo(np.float64).tiny)

            def q_of_z(z):
                """
                Evaluate the Lambert-W transformed variable at axial position z.
                """

                tau_z = 4.0 * kt * z / (U * d)

                # Same beta_tau definition as analytical_efficiency()
                beta_tau = beta - tau_z / b - 1.0

                # Lambert-W argument magnitude:
                # argument = sign * exp(beta_tau) / b
                log_argument_abs = beta_tau - np.log(b)

                if sign > 0.0:
                    # Normal extraction: principal real branch W_0
                    if log_argument_abs < log_max:
                        argument = np.exp(log_argument_abs)
                        q = lambertw(argument, k=0).real
                    else:
                        # Large-positive-argument approximation
                        q = log_argument_abs - np.log(log_argument_abs)

                else:
                    # Inverse permeation: principal real branch W_0
                    if log_argument_abs < log_tiny:
                        argument = 0.0
                    else:
                        argument = -np.exp(log_argument_abs)

                    # Protect against small round-off excursions below -1/e
                    argument = np.clip(argument, -1.0 / np.e, 0.0)
                    q = lambertw(argument, k=0).real

                return float(q)

            def c_bulk(z):
                """
                Molten-salt bulk concentration c_b,l(z).

                This is the local form of the manuscript's cbz_ms_pav equation.
                """

                q = q_of_z(z)
                s = b * (1.0 + q)

                c_b = alpha / 4.0 * (s**2 - 1.0 - 4.0 * Pi_ext)

                # Remove tiny negative values caused only by floating-point error
                if c_b < 0.0 and abs(c_b) < 1.0e-12 * max(c_in, 1.0):
                    c_b = 0.0

                return float(c_b)

            # Integrate the bulk concentration along one pipe
            axial_integral, _ = integrate.quad(
                c_bulk,
                0.0,
                L,
                epsabs=1.0e-12,
                epsrel=1.0e-9,
                limit=200,
            )

        f_H2_to_H = 2.0

        inventory_one_pipe = f_H2_to_H * area_fluid * axial_integral

        self.fluid.inv = inventory_one_pipe * n_pipes
        return self.fluid.inv

    def get_inventory(self, flag_an: bool = True, p_out: float = 0) -> None:
        self.get_solid_inventory(flag_an=flag_an, p_out=p_out)
        self.get_fluid_inventory(flag_an=flag_an, p_out=p_out)
        self.inv = self.fluid.inv + self.membrane.inv
        return
