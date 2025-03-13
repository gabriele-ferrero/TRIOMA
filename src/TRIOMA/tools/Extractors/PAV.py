from TRIOMA.tools.Extractors.PipeSubclasses import (
    Geometry,
    Fluid,
    Membrane,
    FluidMaterial,
    SolidMaterial,
    Turbulator,
    WireCoil,
    CustomTurbulator,
)
from TRIOMA.tools.TriomaClass import TriomaClass
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import lambertw
from scipy import integrate
from typing import Union
from scipy import integrate
from TRIOMA.tools import correlations as corr
import TRIOMA.tools.molten_salts as MS
import TRIOMA.tools.liquid_metals as LM


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
        p_out: float = 1e-15,
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
        self.n_pipes = (self.geometry.n_pipes,)
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

    def custom_update_attribute(self, attr_name: str, new_value: float):
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
        raise ValueError(
            f"'{attr_name}' is not an attribute of {self.__class__.__name__}"
        )

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

    def update_T_prop(self):
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

    def estimate_cost(self, metal_cost=0, fluid_cost=0):
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
        self.pumping_power = (
            self.delta_p * self.get_pipe_flowrate() * self.geometry.n_pipes
        )
        return self.pumping_power

    # def connect_to_component(
    #     self, component2: Union["Component", "BreedingBlanket"] = None
    # ):
    #     """sets the inlet conc of the object component equal to the outlet of self"""
    #     component2.update_attribute("c_in", self.c_out)
    def connect_to_component(self):
        return  ## empty method defined in component_tools.py

    def plot_component(self):
        r_tot = (self.geometry.D) / 2 + self.geometry.thick
        # Create a figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

        # First subplot: two overlapping circles
        circle1 = plt.Circle((0, 0), r_tot, color="#C0C0C0", label="Membrane")
        circle2 = plt.Circle(
            (0, 0), self.geometry.D / 2, color="#87CEEB", label="Fluid"
        )
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
        ax2.arrow(
            0.0, 0.5, 0.1, 0, head_width=0.05, head_length=0.1, fc="black", ec="black"
        )
        # Arrow pointing out of the right side of the rectangle
        ax2.arrow(
            0.8, 0.5, 0.1, 0, head_width=0.05, head_length=0.1, fc="black", ec="black"
        )
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
        Calculates the concentration of the component at the outlet.

        Returns:
            float: The concentration of the component at the outlet.
        """
        if self.fluid.recirculation == 0:
            self.c_out = self.c_in * (1 - self.eff)
        elif self.fluid.recirculation > 0:
            err = 1
            tol = 1e-6
            if self.c_in == 0:
                RaiseError("The inlet concentration is zero")
            c0 = self.c_in
            c_in = c0
            while err > tol:

                c_in1 = c_in
                self.update_attribute("c_in", c_in)
                self.analytical_efficiency()
                self.update_attribute("eff", self.eff_an)
                self.c_out = self.c_in * (1 - self.eff)
                c_in = (self.c_out * self.fluid.recirculation + c0) / (
                    self.fluid.recirculation + 1
                )
                err = abs((c_in - c_in1) / c_in)
        elif self.fluid.recirculation < 0:
            if self.fluid.recirculation <= -1:
                RaiseError(
                    "Bypass(negative recirculation) not valid: it is more than the flowrate"
                )
            if self.c_in == 0:
                RaiseError("The inlet concentration is zero")
            self.c_out = self.c_in * (1 - self.eff) * (
                1 + self.fluid.recirculation
            ) + self.c_in * (-self.fluid.recirculation)
        else:
            RaiseError("Recirculation factor not valid")
        return self.c_out

    def converge_split_HX(
        self,
        tol: float = 1e-3,
        T_in_hot: float = None,
        T_out_hot: float = None,
        T_in_cold: float = None,
        T_out_cold: float = None,
        R_sec: float = None,
        Q: float = None,
        plotvar: bool = False,
        savevar: bool = False,
    ) -> "Circuit":
        """
        Splits the component into N components to better discretize Temperature effects
        Tries to find the optimal number of components to split the component into

        """
        import copy
        from TRIOMA.tools.Circuit import Circuit

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
            axs[1].set_ylabel(
                r"Relative error (%) in $\eta$ with respect to 100 components"
            )
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
        T_in_hot: int = None,
        T_out_hot: int = None,
        T_in_cold: int = None,
        T_out_cold: int = None,
        R_sec: int = 0,
        Q: int = None,
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
            R_cond = np.log(
                (self.fluid.d_Hyd + self.membrane.thick) / self.fluid.d_Hyd
            ) / (2 * np.pi * self.membrane.k)
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
            axes[0].legend(
                ["Primary fluid", "Secondary fluid", "Membrane"], frameon=False
            )
            axes[0].set_ylabel("Temperature [K]")
            axes[0].set_xlabel("Component number")
            axes[0].spines["top"].set_visible(False)
            axes[0].spines["right"].set_visible(False)

            # Second subplot
            axes[1].plot(position_vec, T_vec_p)
            axes[1].plot(position_vec, T_vec_s)
            axes[1].plot(position_vec[:-1], T_vec_membrane)
            axes[1].legend(
                ["Primary fluid", "Secondary fluid", "Membrane"], frameon=False
            )
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

    def get_regime(self, print_var: bool = False):
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

    def get_pipe_flowrate(self):
        """
        Calculates the volumetric flow rate of the component [m^3/s].

        Returns:
            float: The flow rate of the component.
        """
        self.pipe_flowrate = self.fluid.U0 * np.pi * self.fluid.d_Hyd**2 / 4
        return self.pipe_flowrate

    def get_total_flowrate(self):
        """
        Calculates the total flow rate of the component.
        """
        self.get_pipe_flowrate()
        self.flowrate = self.pipe_flowrate * self.geometry.n_pipes
        return self.flowrate * self.geometry.n_pipes

    def define_component_volumes(self):
        """
        Calculates the volumes of the component.
        """
        self.fluid.V = self.geometry.get_fluid_volume()
        self.membrane.V = self.geometry.get_solid_volume()
        self.V = self.fluid.V + self.membrane.V

    def get_adimensionals(self):
        """
        Calculates the adimensional parameters H and W.

        Updates the H and W attributes of the Component object.
        """
        if self.fluid is None:
            print("No fluid selected")
            return
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        match self.fluid.MS:
            case True:
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

    def use_analytical_efficiency(self, p_out=1e-15):
        """Evaluates the analytical efficiency and substitutes it in the efficiency attribute of the component.

        Args:
            L (float): the length of the pipe component
        Returns:
            None
        """
        self.analytical_efficiency(p_out=p_out)
        self.eff = self.eff_an

    def get_efficiency(self, plotvar: bool = False, c_guess: float = None, p_out=1e-15):
        """
        Calculates the efficiency of the component.
        """

        if self.c_in == 0:
            self.c_out = 0
            self.eff = 0
            return

        L_vec = np.linspace(0, self.geometry.L, 100)
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
                    c_guess = self.get_flux(
                        c_vec[i], c_guess=float(self.c_in), p_out=p_out
                    )
            else:
                c_vec[i] = c_vec[
                    i - 1
                ] + f_H2 * self.J_perm * self.fluid.d_Hyd * np.pi * dl**2 / self.fluid.U0 / (
                    np.pi * self.fluid.d_Hyd**2 / 4 * dl
                )
                if isinstance(c_guess, float):
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess, p_out=p_out)
                else:
                    c_guess = self.get_flux(
                        c_vec[i], c_guess=float(self.c_in), p_out=p_out
                    )
        if plotvar:
            plt.plot(L_vec, c_vec)
        self.eff = (self.c_in - c_vec[-1]) / self.c_in

    def analytical_efficiency(self, p_out=1e-15):
        """
        Calculate the analytical efficiency of a component.

        Parameters:
        - p_out: Pressure of H isotope at the outlet of the component. Defaults to 1E-15.

        Returns:
        - eff_an: Analytical efficiency of the component (from Humrickhouse papers)

        This function calculates the analytical efficiency of a component based on the given length (L) of the component.
        It uses various properties of the component, such as membrane properties, fluid properties, and adimensionals.



        If the fluid is a molten salt (MS=True), the analytical efficiency is calculated using the following formula:
        eff_an = 1 - xi * (lambertw(z=np.exp(beta - tau - 1), tol=1e-10) ** 2 + 2 * lambertw(z=np.exp(beta - tau - 1), tol=1e-10)).

        If the fluid is a liquid metal (MS= False), the analytical efficiency is calculated using the following formula:
        eff_an = 1 - np.exp(-tau * zeta / (1 + zeta))* (1-p_in/p_out)^0.5.

        The output of the function is the analytical efficiency of the component as Component.eff_an.
        """
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        self.tau = (
            4 * self.fluid.k_t * self.geometry.L / (self.fluid.U0 * self.fluid.d_Hyd)
        )
        match self.fluid.MS:
            case True:
                self.alpha = (
                    1
                    / (self.fluid.Solubility)
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
                self.xi = self.alpha / self.c_in
                p_in = self.c_in / self.fluid.Solubility
                match (self.xi, self.tau):
                    case (self.xi, self.tau) if self.xi > 1e5:
                        corr_p = 1 - (p_out / p_in)
                        self.eff_an = (1 - np.exp(-self.tau)) * corr_p
                    case (
                        self.xi,
                        self.tau,
                    ) if (
                        self.xi**0.5 < 1e-2 and self.tau > 1 / self.xi**0.5
                    ):
                        corr_p = 1 - (p_out / p_in) ** 0.5
                        self.eff_an = (1 - (1 - self.tau * self.xi**0.5) ** 2) * corr_p
                    case _:
                        e = (self.alpha * p_out * self.fluid.Solubility) ** 0.5
                        f = e / self.alpha
                        delta = (1 / self.xi + 1 + 2 * f) ** 0.5
                        beta = delta + (1 + f) * np.log(delta - 1 - f)
                        max_exp = np.log(np.finfo(np.float64).max)
                        beta_tau = beta - self.tau - 1
                        if beta_tau > max_exp or p_out > 1e-5:
                            # we can use the approximation w=beta_tau-np.log(beta_tau)for the lambert W function but it leads to error up to 40 % in very niche scenarios.

                            def eq(var):
                                cl = var
                                self.alpha = self.xi * self.c_in

                                left = (cl / self.alpha + 1 + 2 * f) ** 0.5 + (
                                    1 + f
                                ) * np.log(
                                    -f + ((cl / self.alpha + 1 + 2 * f) ** 0.5 - 1)
                                )

                                right = beta - self.tau

                                return abs(left - right)

                            p_in = self.c_in / self.fluid.Solubility
                            if (
                                abs(self.p_out * self.fluid.Solubility - self.c_in)
                                / self.c_in
                                < 1e-2
                            ):
                                self.eff_an = 1e-6
                                return
                            lower_bound = min(
                                self.p_out * self.fluid.Solubility, self.c_in
                            )
                            upper_bound = max(
                                self.p_out * self.fluid.Solubility, self.c_in
                            )
                            cl = minimize(
                                eq,
                                x0=(lower_bound + upper_bound) / 2,
                                method="Powell",
                                bounds=[(lower_bound, upper_bound)],
                                tol=1e-7,
                            ).x[0]
                            # corr_p=1-(p_out/p_in)
                            self.eff_an = 1 - (cl / self.c_in)
                            return
                        else:
                            z = np.exp(beta_tau)
                            w = lambertw(z, tol=1e-10)
                            self.eff_an = 1 - self.xi * (w**2 + 2 * w)
                            if self.eff_an.imag != 0:
                                raise ValueError(
                                    "self.eff_an has a non-zero imaginary part"
                                )
                            else:
                                self.eff_an = self.eff_an.real  # get rid of 0*j
            case False:  # Liquid Metal
                self.zeta = (2 * self.membrane.K_S * self.membrane.D) / (
                    self.fluid.k_t
                    * self.fluid.Solubility
                    * self.fluid.d_Hyd
                    * np.log(
                        (self.fluid.d_Hyd + 2 * self.membrane.thick) / self.fluid.d_Hyd
                    )
                )
                p_in = (self.c_in / self.fluid.Solubility) ** 2
                corr_p = 1 - (p_out / p_in) ** 0.5

                self.eff_an = (
                    1 - np.exp(-self.tau * self.zeta / (1 + self.zeta))
                ) * corr_p

    def get_flux(self, c: float = None, c_guess: float = 1e-9, p_out=1e-15):
        """
        Calculates the Tritium flux of the component.
        It can make some approximations based on W and H to make the solver faster

        Args:
            c (float): The concentration in the component fluid.

        Returns:
            float: The permeation flux.

        """
        if not isinstance(c, float):
            print(c)
            raise ValueError("Input 'c' must be a non-empty numpy array")

        if not isinstance(c_guess, float):
            raise ValueError("c_guess must be a float")
        self.get_adimensionals()
        if self.fluid.MS:
            if self.W > 10:
                # DIFFUSION LIMITED V // Surface limited X
                if self.H / self.W > 1000:
                    # Mass transport limited V // Diffusion limited X
                    self.J_perm = (
                        -2 * self.fluid.k_t * (c - p_out * self.fluid.Solubility)
                    )  ## MS factor
                elif self.H / self.W < 0.0001:
                    # Diffusion limited V // Mass Transport limited X
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
                        * ((c / self.fluid.Solubility) ** 0.5 - p_out**0.5)
                    )
                else:
                    # Mixed regime mass transport diffusion
                    def equations(vars):
                        if vars.size == 0:
                            return upper_bound
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
                                * ((c_wl / self.fluid.Solubility) ** 0.5 - p_out**0.5)
                            )
                        )
                        return abs(J_diff - J_mt)

                    if isinstance(c_guess, float):

                        initial_guess = c_guess
                    else:
                        ValueError("c_guess must be a float")
                    initial_guess = c_guess
                    min_upper_bound = 1e-4  # Set a minimum value for the upper bound
                    upper_bound = max(c * (1 + 1e-4), min_upper_bound)
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, upper_bound),
                        ],
                        tol=1e-8,
                        options={
                            "maxiter": int(1e7),
                        },
                    )
                    self.J_perm = (
                        -2 * self.fluid.k_t * (c - solution.x[0])
                    )  ## MS factor
                    return float(solution.x[0])
            elif self.W < 0.1:
                # Surface limited V // Diffusion Limited X
                if self.H > 100:
                    # Mass transport limited V // Surface limited X
                    self.J_perm = (
                        -2 * self.fluid.k_t * (c - p_out * self.fluid.Solubility)
                    )  ## MS factor
                elif self.H < 0.01:
                    # Surface limited V // Mass Transport limited X
                    self.J_perm = -self.membrane.k_d * (c / self.fluid.Solubility)
                else:
                    # Mixed regime mass transfer surface
                    def equations(vars):
                        if vars.size == 0:
                            return upper_bound
                        c_wl = vars
                        c_bl = c
                        J_mt = 2 * self.fluid.k_t * (c_bl - c_wl)  ## MS factor
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )

                        return abs(J_mt - J_surf)

                    initial_guess = [c * 1e-1]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c * (1 + 1e-4)),
                        ],
                        tol=1e-8,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_bl = c
                    c_wl = solution.x[0]
                    self.J_perm = 2 * self.fluid.k_t * (c_bl - c_wl)  ## MS factor
                    return float(solution.x[0])
            else:
                # Mixed Diffusion Surface
                if self.H / self.W > 1000:
                    # Mass transport limited V // Mixed Surface Diffusion Limited X
                    self.J_perm = (
                        -2 * self.fluid.k_t * (c - p_out * self.fluid.Solubility)
                    )  ## MS factor
                elif self.H / self.W < 0.0001:
                    # Mixed Diffusion Surface
                    def equations(vars):
                        if vars.size == 0:
                            return upper_bound
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
                                * ((c_wl / self.fluid.Solubility) ** 0.5 - p_out**0.5)
                            )
                        )

                        return abs(J_diff - J_surf)

                    if c_guess is None:
                        initial_guess = [c / 2]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (1e-14, c),
                        ],
                        tol=1e-7,
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
                        * (
                            self.membrane.K_S * (c_wall / self.fluid.Solubility) ** 0.5
                            - p_out**0.5
                        )
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
                        ) - self.membrane.k_d * self.membrane.K_S**2 * (c_ws**2)
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
                            * ((self.membrane.K_S * c_ws) - p_out**0.5)
                        )
                        eq1 = abs(J_mt - J_d)
                        eq2 = abs(J_mt - J_diff)
                        eq3 = abs(J_d - J_diff)

                        return eq1 + eq2 + eq3

                    initial_guess = [(2 * c / 3), (c / 3)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[(0, c), (0, c)],
                        tol=1e-8,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_wl = solution.x[0]
                    self.J_perm = 2 * self.fluid.k_t * (c - c_wl)
                    result = float(solution.x[0])
                    if np.isscalar(result):
                        return result

        else:
            if self.W > 10:
                # DIFFUSION LIMITED V // Surface limited X
                if self.H / self.W > 1000:
                    # Mass transport limited V // Diffusion limited X
                    self.J_perm = -self.fluid.k_t * (
                        c - p_out**0.5 * self.fluid.Solubility
                    )  ## LM factor
                elif self.H / self.W < 0.0001:
                    # Diffusion limited V // Mass Transport limited X
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
                        * (self.membrane.K_S * (c / self.fluid.Solubility - p_out**0.5))
                    )
                else:
                    # Mixed regime mass transport diffusion
                    def equations(vars):
                        c_wl = vars
                        J_mt = self.fluid.k_t * (c - c_wl)  ## LM factor
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
                                * (c_wl / self.fluid.Solubility - p_out**0.5)
                            )
                        )
                        return abs((J_diff - J_mt))

                    if c_guess is None:
                        initial_guess = [c * 1e-2]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c * (1 + 1e-4)),
                        ],
                        tol=1e-8,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    self.J_perm = -self.fluid.k_t * (c - solution.x[0])  ## LM factor
                    return float(solution.x[0])
            elif self.W < 0.1:
                # Surface limited V // Diffusion Limited X
                if self.H > 100:
                    # Mass transport limited V // Surface limited X
                    self.J_perm = -self.fluid.k_t * (
                        c - p_out**0.5 * self.fluid.Solubility
                    )  ## LM factor
                elif self.H < 0.01:
                    # Surface limited V // Mass Transport limited X
                    self.J_perm = -self.membrane.k_d * (c / self.fluid.Solubility)
                else:
                    # Mixed regime mass transfer surface
                    def equations(vars):
                        c_wl = vars
                        c_bl = c
                        J_mt = self.fluid.k_t * (
                            c_bl - c_wl - p_out**0.5 * self.fluid.Solubility
                        )  ## LM factor
                        J_surf = (
                            self.membrane.k_d * (c_bl / self.fluid.Solubility)
                            - self.membrane.k_d * self.membrane.K_S**2 * c_wl**2
                        )

                        return abs(J_mt - J_surf)

                    initial_guess = [c * 1e-1]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c * (1 + 1e-4)),
                        ],
                        tol=1e-8,
                        options={
                            "maxiter": int(1e6),
                        },
                    )
                    c_bl = c
                    c_wl = solution.x[0]
                    self.J_perm = self.fluid.k_t * (
                        c_bl - c_wl - p_out * self.fluid.Solubility
                    )  ## LM factor
                    return float(solution.x[0])
            else:
                if self.H / self.W < 0.0001:
                    # Mass transport limited V // Mixed Surface Diffusion  X
                    self.J_perm = -self.fluid.k_t * (
                        c - p_out**0.5 * self.fluid.Solubility
                    )  ## LM factor
                elif self.H / self.W > 1000:
                    # Mixed Diffusion Surface V // Mass Transport limited X
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
                                * (c_wl / self.fluid.Solubility - p_out**0.5)
                            )
                        )

                        return abs(J_diff - J_surf)

                    if c_guess is None:
                        initial_guess = [c / 2]
                    else:
                        initial_guess = c_guess
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (1e-14, c),
                        ],
                        tol=1e-8,
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
                        * (
                            self.membrane.K_S
                            * (c_wall / self.fluid.Solubility - p_out**0.5)
                        )
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
                            * ((self.membrane.K_S * c_ws) - p_out**0.5)
                        )
                        eq1 = abs(J_mt - J_d)
                        eq2 = abs(J_mt - J_diff)
                        eq3 = abs(J_d - J_diff)

                        return eq1 + eq2 + eq3

                    initial_guess = [(2 * c / 3), (c / 3)]
                    solution = minimize(
                        equations,
                        initial_guess,
                        method="Powell",
                        bounds=[
                            (0, c * (1 + 1e-4)),
                        ],
                        tol=1e-8,
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
        if self.geometry.turbulator is None:
            h_prim = corr.get_h_from_Nu(
                corr.Nu_DittusBoelter(Re, Pr), self.fluid.k, self.fluid.d_Hyd
            )
        else:
            match self.geometry.turbulator.turbulator_type:
                case "TwistedTape":
                    print(
                        str(self.geometry.turbulator.turbulator_type)
                        + " is not implemented yet"
                    )
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

    def analytical_solid_inventory(self, p_out: float = 0):
        if self.fluid.k_t is None:

            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        match self.fluid.MS:
            case False:

                def integralfun(r):
                    return (
                        1
                        / 4
                        * r**2
                        * (
                            2 * np.log(r / (self.geometry.D / 2 + self.geometry.thick))
                            - 1
                        )
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
                        * np.log(
                            (self.fluid.d_Hyd + 2 * self.membrane.thick)
                            / self.fluid.d_Hyd
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
                            (self.fluid.d_Hyd + 2 * self.membrane.thick)
                            / self.fluid.d_Hyd
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
                K = (
                    -2
                    * np.pi
                    * (
                        self.c_in
                        / (dimless2 / self.fluid.k_t + 1)
                        / self.fluid.Solubility
                        * self.membrane.K_S
                    )
                    / np.log(
                        (self.geometry.D / 2 + self.geometry.thick)
                        / (self.geometry.D / 2)
                    )
                )
                L_factor = (np.exp(L_ch * self.geometry.L) - 1) / L_ch
                K = K * L_factor
                integral = (
                    K * integralfun(self.geometry.D / 2 + self.geometry.thick)
                    - K * integralfun(self.geometry.D / 2)
                ) + self.geometry.L * p_out**0.5 * self.membrane.K_S * (
                    circle(self.geometry.D / 2 + self.geometry.thick)
                    - circle(self.geometry.D / 2)
                )
                inventory = integral
                self.membrane.inv = inventory
                return inventory
            case True:

                def ms_integral(self, p_out: float = 0, L: float = 0):
                    if self.tau is None or self.xi is None or self.alpha is None:
                        self.analytical_efficiency(p_out=p_out)
                    beta = (1 / self.xi + 1) ** 0.5 + np.log(
                        (1 / self.xi + 1) ** 0.5 - 1
                    )
                    max_exp = np.log(np.finfo(np.float64).max)
                    beta_tau = beta - self.tau - 1
                    if beta_tau > max_exp:
                        w = beta_tau - np.log(beta_tau)
                        w2 = -beta_tau

                    else:
                        z = np.exp(beta_tau)
                        z2 = np.exp(-beta_tau)
                        w = lambertw(z, tol=1e-10)
                        w2 = lambertw(z2, tol=1e-10)
                        if w.imag != 0:
                            raise ValueError(
                                "self.eff_an has a non-zero imaginary part"
                            )
                        if w2.imag != 0:
                            raise ValueError(
                                "self.eff_an has a non-zero imaginary part"
                            )
                        w = w.real
                        w2 = w2.real
                    c_ext = p_out**0.5 * self.membrane.K_S
                    conv = (
                        self.c_in / self.fluid.Solubility
                    ) ** 0.5 * self.membrane.K_S
                    c_w_l = self.alpha * (w**2 + 2 * w) + self.alpha * (
                        2 - 2 * ((w**2 + 2 * w) + 1) ** 0.5
                    )
                    K = (
                        self.alpha**0.5
                        / self.fluid.Solubility**0.5
                        * (
                            -beta_tau
                            * (w**2 - w + 1)
                            / (
                                4
                                * self.fluid.k_t
                                / (self.fluid.U0 * self.fluid.d_Hyd)
                                * w2
                            )
                        )
                        * self.membrane.K_S
                    )

                    def integralfun(r):
                        return (
                            1
                            / 4
                            * r**2
                            * (
                                2
                                * np.log(
                                    r / (self.geometry.D / 2 + self.geometry.thick)
                                )
                                - 1
                            )
                        )

                    integral = K * integralfun(
                        self.geometry.D / 2 + self.geometry.thick
                    ) - K * integralfun(self.geometry.D / 2)
                    return integral

                def p_out_term(self, p_out):
                    def circle(r):
                        return np.pi * r**2

                    add = (
                        self.geometry.L
                        * p_out**0.5
                        * self.membrane.K_S
                        * (
                            circle(self.geometry.D / 2 + self.geometry.thick)
                            - circle(self.geometry.D / 2)
                        )
                    )

                    return add

                inv = (
                    ms_integral(self=self, L=self.geometry.L, p_out=p_out)
                    - ms_integral(self=self, L=0, p_out=p_out)
                    + p_out_term(self, p_out)
                )
                self.membrane.inv = inv * self.geometry.n_pipes
                return inv

    def get_solid_inventory(self, p_out: float = 0, flag_an=False):
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
                                (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                / self.fluid.d_Hyd
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
                                (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                / self.fluid.d_Hyd
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
                        (
                            -(c_w - c_ext) * np.log(r / r_out) / np.log(r_out / r_in)
                            + c_ext
                        )
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
                                    (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                    / self.fluid.d_Hyd
                                )
                            )
                        )
                        ** 2
                    )

                    beta = (1 / self.xi + 1) ** 0.5 + np.log(
                        (1 / self.xi + 1) ** 0.5 - 1
                    )
                    max_exp = np.log(np.finfo(np.float64).max)
                    beta_tau = beta - tau - 1
                    if beta_tau > max_exp:

                        w = beta_tau - np.log(beta_tau)

                    else:
                        z = np.exp(beta_tau)
                        w = lambertw(z, tol=1e-10)
                        if w.imag != 0:
                            raise ValueError(
                                "self.eff_an has a non-zero imaginary part"
                            )
                        w = w.real
                    alpha = (
                        1
                        / self.fluid.Solubility
                        * (
                            (
                                0.5  ## TODO: Check this
                                * self.membrane.D
                                * self.membrane.K_S
                            )
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
                    c_ext = p_out**0.5 * self.membrane.K_S
                    conv = (
                        self.c_in / self.fluid.Solubility
                    ) ** 0.5 * self.membrane.K_S
                    c_w_l = (
                        alpha * (w**2 + 2 * w)
                        + alpha
                        * (2 - 2 * ((w**2 + 2 * w) + 1) ** 0.5)  ## TODO: Check this
                        + c_ext
                    )

                    if c_w_l < 0:
                        c_w_l = 1e-17
                    return (
                        (
                            -np.log(r / r_out)
                            / np.log(r_out / r_in)
                            * (
                                (alpha / self.fluid.Solubility) ** 0.5
                                * w
                                * self.membrane.K_S
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

    def analytical_fluid_inventory(self, p_out=0):
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
                        * np.log(
                            (self.fluid.d_Hyd + 2 * self.membrane.thick)
                            / self.fluid.d_Hyd
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
                            (self.fluid.d_Hyd + 2 * self.membrane.thick)
                            / self.fluid.d_Hyd
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
                c_ext = p_out**0.5 * self.fluid.Solubility
                K = self.c_in - c_ext

                L_factor = (np.exp(L_ch * self.geometry.L) - 1) / L_ch
                K = K * L_factor
                integral = (K) * circle(self.geometry.D / 2) + c_ext * circle(
                    self.geometry.D / 2
                ) * self.geometry.L
                inventory = integral
                self.fluid.inv = inventory * self.geometry.n_pipes
                return inventory
            case True:
                print("MS fluid integration is done numerically")
                self.get_fluid_inventory(flag_an=False, p_out=p_out)

    def get_fluid_inventory(self, flag_an=False, p_out=0):
        if flag_an == True:
            return self.analytical_fluid_inventory(p_out=p_out)
        r_in = self.fluid.d_Hyd / 2

        L_min = 0
        L_max = self.geometry.L
        N = 100

        def integrand(L):
            if self.fluid.k_t is None:
                self.fluid.get_kt(turbulator=self.geometry.turbulator)
            match self.fluid.MS:
                case False:
                    dimless = (
                        2
                        * self.membrane.D
                        * self.membrane.K_S
                        / (
                            self.fluid.k_t
                            * self.fluid.Solubility
                            * self.fluid.d_Hyd
                            * np.log(
                                (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                / self.fluid.d_Hyd
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
                    c_ext = p_out**0.5 * self.fluid.Solubility
                    return (self.c_in - c_ext) * np.exp(L_ch * L) + c_ext
                case True:
                    if self.tau is None or self.xi is None:
                        self.analytical_efficiency(p_out=p_out)
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
                                    (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                    / self.fluid.d_Hyd
                                )
                            )
                        )
                        ** 2
                    )

                    beta = (1 / self.xi + 1) ** 0.5 + np.log(
                        (1 / self.xi + 1) ** 0.5 - 1
                    )
                    max_exp = np.log(np.finfo(np.float64).max)
                    beta_tau = beta - tau - 1
                    if beta_tau > max_exp:
                        # print(
                        #     "Warning: Overflow encountered in exp, input too large.Approximation triggered"
                        # )

                        w = beta_tau - np.log(beta_tau)
                    else:
                        z = np.exp(beta_tau)
                        w = lambertw(z, tol=1e-10)
                        if w.imag != 0:
                            raise ValueError(
                                "self.eff_an has a non-zero imaginary part"
                            )
                        w = w.real
                    alpha = (
                        1
                        / self.fluid.Solubility
                        * (
                            (
                                0.5 * self.membrane.D * self.membrane.K_S
                            )  ## TODO: Check this
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
                    conv = (
                        self.c_in / self.fluid.Solubility
                    ) ** 0.5 * self.membrane.K_S
                    c_w_l = alpha * (w**2 + 2 * w)
                    return c_w_l

        result, err = integrate.nquad(integrand, [[L_min, L_max]])
        self.fluid.inv = result * np.pi * r_in**2 * self.geometry.n_pipes
        return self.fluid.inv

    def get_inventory(self, flag_an=True, p_out=0):
        self.get_solid_inventory(flag_an=flag_an, p_out=p_out)
        self.get_fluid_inventory(flag_an=flag_an, p_out=p_out)
        self.inv = self.fluid.inv + self.membrane.inv
        return
