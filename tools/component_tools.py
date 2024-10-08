from ast import Raise
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
from typing import Union
import matplotlib.pyplot as plt
from scipy import integrate


def print_class_variables(instance, variable_names=None, tab: int = 0):
    """
    Prints specified variables of a class. If a variable is a class itself,
    calls this function recursively for the internal class. If variable_names is None,
    prints all variables.

    Args:
        instance: The class instance.
        variable_names (list of str, optional): Names of the variables to print. Prints all if None.
    """
    built_in_types = [
        Component,
        Fluid,
        Membrane,
        FluidMaterial,
        SolidMaterial,
        BreedingBlanket,
        Geometry,
        Turbulator,
    ]
    indent = "    " * tab  # Define the indentation as four spaces per tab level
    for attr_name, attr_value in instance.__dict__.items():
        if variable_names is None or attr_name.lower() == variable_names.lower():
            if type(attr_value) in built_in_types:
                tab += 1
                print(
                    f"{indent}{attr_name} is a {type(attr_value)} class, printing its variables:"
                )
                print_class_variables(attr_value, variable_names, tab=tab)
                tab -= 1
            else:
                print(f"{indent}{attr_name}: {attr_value}")


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
        for attr, value in instance.__dict__.items():
            if isinstance(value, object) and hasattr(value, attr_name):
                setattr(value, attr_name, new_value)
                return
        raise ValueError(
            f"'{attr_name}' is not an attribute of {instance.__class__.__name__}"
        )


class Geometry:
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
        n_pipes: int = 1,
        turbulator: Union["Turbulator"] = None,
    ):
        self.L = L
        self.D = D
        self.thick = thick
        self.n_pipes = n_pipes
        self.turbulator = turbulator

    def update_attribute(self, attr_name: str, new_value: float):
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


class Circuit:
    """
    Represent a circuit of components connected in series

    This class represents a circuit consisting of multiple components connected in series. It provides methods to update attributes, add components, calculate circuit efficiency, and plot the circuit.

    Attributes:
        components (list): A list of components in the circuit.
        closed (bool): A boolean indicating whether the circuit is a closed loop or not.

    Methods:
        update_attribute(attr_name, new_value):
            Updates the value of the specified attribute.

        add_component(component):
            Adds a component to the circuit.

        get_eff_circuit():
            Calculates the efficiency of the circuit.

        get_gains_and_losses():
            Calculates the gains and losses of the circuit.

        plot_circuit():
            Plots the circuit using matplotlib.

    Example usage:
        circuit = Circuit()
        circuit.add_component(component1)
        circuit.add_component(component2)
        circuit.get_eff_circuit()
        circuit.plot_circuit()
    """

    def __init__(
        self,
        components: list = None,
        closed: bool = False,
    ):
        vec_components = []
        if components is not None:
            for element in components:
                if isinstance(element, Union[Component, BreedingBlanket]):
                    vec_components.append(element)
                elif isinstance(element, Circuit):
                    for comp in element.components:
                        vec_components.append(comp)
                else:
                    raise ValueError("Invalid component type")
        self.components = vec_components
        self.closed = closed

    def update_attribute(self, attr_name: str, new_value: Union[float, bool]):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def add_component(
        self, component: Union["Component", "BreedingBlanket", "Circuit"]
    ):
        """
        Adds a component to the circuit.

        Args:
            component (Component): The component to add.
        """
        if isinstance(component, Circuit):
            for comp in component.components:
                self.components.append(comp)
        else:
            self.components.append(component)

    def get_eff_circuit(self):
        """
        Calculates the efficiency of the circuit based on the components present.

        Returns:
            circtuit.eff (float): The efficiency of the circuit.

        Raises:
            None

        Example Usage:
            circuit.get_eff_circuit()

        """
        ind = None
        for i, component in enumerate(self.components):
            if isinstance(component, BreedingBlanket):
                ind = i
                break  # Assuming there's only one BreedingBlanket

        if ind is not None:
            # Move the element at index `ind` to the first position
            self.components = (
                [self.components[ind]]
                + self.components[:ind]
                + self.components[ind + 1 :]
            )

        for i, component in enumerate(self.components):
            if isinstance(component, Component):
                component.use_analytical_efficiency()
                component.outlet_c_comp()
            if i != len(self.components) - 1:
                component.connect_to_component(self.components[i + 1])
        eff_circuit = (
            self.components[1].c_in - self.components[-1].c_out
        ) / self.components[1].c_in
        self.eff = eff_circuit

    def get_gains_and_losses(self):
        """
        Calculates the gains and losses of the components in the circuit.

        Returns:
            circuit.extraction_perc (float): The extraction percentage of the circuit.
            circuit.loss_perc (float): The loss percentage of the circuit.

        """

        gains = 0
        losses = 0
        flag_bb = 0
        for i, component in enumerate(self.components):
            diff = component.c_in - component.c_out
            if isinstance(component, Component):
                if component.loss == False:
                    gains += diff
                else:
                    losses += diff
        for i, component in enumerate(self.components):

            if isinstance(component, BreedingBlanket):
                ind = i
                if flag_bb != 0:
                    print("There are more BB!")
                flag_bb = 1
        if ind != 0 and ind != len(self.components):
            eff_circuit = (
                self.components[ind + 1].c_in - self.components[ind - 1].c_out
            ) / self.components[ind + 1].c_in
            self.extraction_perc = gains / self.components[ind + 1].c_in / eff_circuit
            self.loss_perc = losses / self.components[ind + 1].c_in / eff_circuit
        elif ind == 0:
            eff_circuit = (
                self.components[ind + 1].c_in - self.components[-1].c_out
            ) / self.components[ind + 1].c_in
            self.extraction_perc = gains / self.components[ind + 1].c_in / eff_circuit
            self.loss_perc = losses / self.components[ind + 1].c_in / eff_circuit
        elif ind == len(self.components):
            eff_circuit = (
                self.components[0].c_in - self.components[ind - 1].c_out
            ) / self.components[0].c_in
            self.extraction_perc = gains / self.components[0].c_in / eff_circuit
            self.loss_perc = losses / self.components[0].c_in / eff_circuit
        self.eff = eff_circuit

    def plot_circuit(self):
        """
        Plot the circuit diagram for the components in the circuit.

        This function uses matplotlib to create a circuit diagram for the components in the circuit.
        Each component is represented by a rectangle with arrows indicating the flow direction.
        The color of the rectangle represents the position of the component in the circuit.

        Returns:
            None
        """
        from matplotlib.colors import LinearSegmentedColormap

        # Define the red-to-blue colormap
        red_to_blue = LinearSegmentedColormap.from_list("RedToBlue", ["red", "blue"])
        num_components = len(self.components)
        if num_components < 10:
            num_rows = 1
            num_columns = num_components
        else:
            num_rows = (num_components // 10) + (1 if num_components % 10 != 0 else 0)
            num_columns = 10

        fig, axs = plt.subplots(
            num_rows + 1, num_columns, figsize=(5 * num_columns, 4 * num_rows)
        )
        fig.subplots_adjust(hspace=0)
        # for component in self.components:
        #     component.plot_component()
        for i, component in enumerate(self.components):
            color = red_to_blue(i / num_components)
            if isinstance(component, Component):
                rectangle = plt.Rectangle(
                    (0.2, 0.3), 0.55, 0.4, edgecolor="black", facecolor=color, alpha=0.5
                )
                axs[i // num_columns, i % num_columns].add_patch(rectangle)
                # Arrow pointing to the left side of the rectangle
                axs[i // num_columns, i % num_columns].arrow(
                    0.0,
                    0.5,
                    0.1,
                    0,
                    head_width=0.05,
                    head_length=0.1,
                    fc="black",
                    ec="black",
                )
                # Arrow pointing out of the right side of the rectangle
                axs[i // num_columns, i % num_columns].arrow(
                    0.8,
                    0.5,
                    0.1,
                    0,
                    head_width=0.05,
                    head_length=0.1,
                    fc="black",
                    ec="black",
                )
                axs[i // num_columns, i % num_columns].set_aspect("equal")
                axs[i // num_columns, i % num_columns].set_xlim(0, 1)
                axs[i // num_columns, i % num_columns].set_ylim(0, 1)
                axs[i // num_columns, i % num_columns].text(
                    0.15,
                    0.3,
                    f"L={component.geometry.L:.3g} m",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.15,
                    0.4,
                    f"T={component.fluid.T:.6g}K",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.15,
                    0.6,
                    f"$c_{{in}}={component.c_in:.4g}  mol/m^3$",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.5,
                    0.7,
                    f"velocity={component.fluid.U0:.2g} m/s",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.5,
                    0.4,
                    f"eff={component.eff*100:.2g}%",
                    color="black",
                    ha="center",
                    va="center",
                )

                axs[i // num_columns, i % num_columns].text(
                    0.9,
                    0.3,
                    rf"$c_{{out}}={component.c_out:.4g} mol/m^3$",
                    color="black",
                    ha="center",
                    va="center",
                )
            elif isinstance(component, BreedingBlanket):

                rectangle = plt.Rectangle(
                    (0.2, 0.3),
                    0.55,
                    0.4,
                    edgecolor="black",
                    facecolor="green",
                    alpha=0.5,
                )
                axs[i // num_columns, i % num_columns].add_patch(rectangle)
                # Arrow pointing to the left side of the rectangle
                axs[i // num_columns, i % num_columns].arrow(
                    0.5,
                    0.7,
                    0,
                    0.1,
                    head_width=0.05,
                    head_length=0.1,
                    fc="black",
                    ec="black",
                )
                # Arrow pointing out of the right side of the rectangle
                axs[i // num_columns, i % num_columns].arrow(
                    0.5,
                    0.1,
                    0.0,
                    0.1,
                    head_width=0.05,
                    head_length=0.1,
                    fc="black",
                    ec="black",
                )
                axs[i // num_columns, i % num_columns].set_aspect("equal")
                axs[i // num_columns, i % num_columns].set_xlim(0, 1)
                axs[i // num_columns, i % num_columns].set_ylim(0, 1)
                if component.name is None:
                    axs[i // num_columns, i % num_columns].set_title("Component ")
                else:
                    axs[i // num_columns, i % num_columns].set_title(component.name)
                axs[i // num_columns, i % num_columns].text(
                    0.7,
                    0.8,
                    r"$T_o$=" + str(component.T_out) + " K",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.7,
                    0.2,
                    r"$T_i$=" + str(component.T_in) + " K",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.3,
                    0.2,
                    f"$c_i$={component.c_in:.4g} $mol/m^3$",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.5,
                    0.6,
                    f"Q={component.Q/1E6:.3g} MW",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.5,
                    0.4,
                    f"TBR={component.TBR:.3g}",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].text(
                    0.3,
                    0.8,
                    rf"$c_o$={component.c_out:.4g}$mol/m^3$",
                    color="black",
                    ha="center",
                    va="center",
                )
                axs[i // num_columns, i % num_columns].axis("off")
                # Display the plot
                fig.tight_layout()
            else:
                raise ValueError("Component not recognized")
        for row in axs:
            for ax in row:
                ax.axis("off")
                ax.set_ylim(0.2, 0.8)
        return fig

    def solve_circuit(self, tol=1e-6):
        """
        Solve the circuit by calculating the concentration of the components at the outlet.
        If the circuit is a closed loop, the concentration of the first component is set to the concentration of the last component until the stationary regime is reached.

        """
        err = 1
        flag = 0
        flag_bb = 0
        for i, component in enumerate(self.components):

            if isinstance(component, BreedingBlanket):
                ind = i
                if flag_bb != 0:
                    print("There are more BB!")
                flag_bb = 1
        while flag == 0:
            for i, component in enumerate(self.components):

                if isinstance(component, Component):
                    component.use_analytical_efficiency()
                    component.outlet_c_comp()
                    if i != len(self.components) - 1:
                        component.connect_to_component(self.components[i + 1])
                if isinstance(component, BreedingBlanket):
                    component.get_cout()
                    if i != len(self.components) - 1:
                        component.connect_to_component(self.components[i + 1])
            if self.closed == True:
                err = (
                    abs(self.components[0].c_in - self.components[-1].c_out)
                    / self.components[0].c_in
                )
                if err < tol:
                    flag = 1
            else:
                flag = 1
            self.components[-1].connect_to_component(self.components[0])

    def inspect_circuit(self, name=None):
        """
        Inspects the circuit components.

        Parameters:
            name (str, optional): The name of the component to inspect. If not provided, all components will be inspected.

        Returns:
            None

        """
        if name is None:
            for component in self.components:
                component.inspect()
        else:
            for component in self.components:
                if component.name == name:
                    component.inspect()


class Component:
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
        n_pipes: int = 1,
        fluid: "Fluid" = None,
        membrane: "Membrane" = None,
        name: str = None,
        loss: bool = False,
        inv: float = None,
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
        self.n_pipes = (n_pipes,)
        self.fluid = fluid
        self.membrane = membrane
        self.name = name
        self.loss = loss
        self.inv = inv
        # if (
        #     isinstance(self.fluid, Fluid)
        #     and isinstance(self.membrane, Membrane)
        #     and isinstance(self.geometry, Geometry)
        # ):
        #     if self.membrane.thick != self.geometry.thick:
        #         print("overwriting Membrane thickness with Geometry thickness")
        #     if self.fluid.d_Hyd != self.geometry.D:
        #         print("overwriting Fluid Hydraulic diameter with Geometry Diameter")
        # self.membrane.thick = self.geometry.thick
        # self.fluid.d_Hyd = self.geometry.D

    def update_attribute(
        self,
        attr_name: str = None,
        new_value: Union[float, "Fluid", "Membrane", "Geometry"] = None,
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

    def connect_to_component(
        self, component2: Union["Component", "BreedingBlanket"] = None
    ):
        """sets the inlet conc of the object component equal to the outlet of self"""
        component2.update_attribute("c_in", self.c_out)

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
        self.c_out = self.c_in * (1 - self.eff)

    def split_HX(
        self,
        N: int = 25,
        T_in_hot: int = None,
        T_out_hot: int = None,
        T_in_cold: int = None,
        T_out_cold: int = None,
        R_sec: int = None,
        Q: int = None,
        plotvar: bool = False,
    ) -> "Circuit":
        """
        Splits the component into N components to better discretize Temperature effects
        """
        import copy

        deltaTML = corr.get_deltaTML(T_in_hot, T_out_hot, T_in_cold, T_out_cold)
        self.get_global_HX_coeff(R_sec)
        L_tot = corr.get_length_HX(
            deltaTML=deltaTML, d_hyd=self.geometry.D, U=self.U, Q=Q
        )
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
                    loss=copy.deepcopy(self.loss),
                )
            )
        T_vec_p = np.linspace(T_in_hot, T_out_hot, N)
        L_vec = []
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

            component.get_global_HX_coeff(R_sec)
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

        if plotvar == True:
            plt.plot(T_vec_p)
            plt.plot(T_vec_s)
            x_values = np.arange(len(T_vec_membrane)) + 0.5
            plt.plot(x_values, T_vec_membrane)
            plt.legend(["Primary fluid", "Secondary fluid", "Membrane"])
            plt.show()

        return circuit

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
    ) -> "Circuit":
        """
        Splits the component into N components to better discretize Temperature effects
        Tries to find the optimal number of components to split the component into

        """
        import copy

        eff_v = []
        for N in range(10, 101, 10):
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
        x_values = range(10, 101, 10)
        if plotvar == True:
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))

            # First subplot
            axs[0].plot(x_values, eff_v)
            axs[0].set_xlabel("Number of components")
            axs[0].set_ylabel("Efficiency")

            # Second subplot
            axs[1].semilogy(x_values, abs(eff_v - eff_v[-1]) / eff_v[-1] * 100)
            axs[1].set_xlabel("Number of components")
            axs[1].set_ylabel(
                f"Relative error in efficiency (%) with respect to {100} components"
            )

            plt.tight_layout()
            plt.show()

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
        if self.fluid is not None:
            if self.fluid.k_t is None:
                self.fluid.get_kt(turbulator=self.geometry.turbulator)
            if self.fluid.MS == True:
                if self.membrane is not None:

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
                else:
                    return "No membrane selected"
            else:
                if self.membrane is not None:
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
                else:
                    return "No membrane selected"
        else:
            return "No fluid selected"

    def get_pipe_flowrate(self):
        """
        Calculates the volumetric flow rate of the component [m^3/s].

        Returns:
            float: The flow rate of the component.
        """
        self.pipe_flowrate = self.fluid.U0 * np.pi * self.fluid.d_Hyd**2 / 4
        return self.fluid.U0 * np.pi * self.fluid.d_Hyd**2 / 4

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
        if self.fluid.k_t is None:
            self.fluid.get_kt(turbulator=self.geometry.turbulator)
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

    def use_analytical_efficiency(self):
        """Evaluates the analytical efficiency and substitutes it in the efficiency attribute of the component.

        Args:
            L (float): the length of the pipe component
        Returns:
            None
        """
        self.analytical_efficiency()
        self.eff = self.eff_an

    def get_efficiency(self, plotvar: bool = False, c_guess: float = None):
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
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=float(self.c_in))
            else:
                c_vec[i] = c_vec[
                    i - 1
                ] + f_H2 * self.J_perm * self.fluid.d_Hyd * np.pi * dl**2 / self.fluid.U0 / (
                    np.pi * self.fluid.d_Hyd**2 / 4 * dl
                )
                if isinstance(c_guess, float):
                    c_guess = self.get_flux(c_vec[i], c_guess=c_guess)
                else:
                    c_guess = self.get_flux(c_vec[i], c_guess=float(self.c_in))
        if plotvar:
            plt.plot(L_vec, c_vec)
        self.eff = (self.c_in - c_vec[-1]) / self.c_in

    def analytical_efficiency(self):
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
        if self.fluid.k_t is None:
            self.fluid.get_kt(turbulator=self.geometry.turbulator)
        self.tau = (
            4 * self.fluid.k_t * self.geometry.L / (self.fluid.U0 * self.fluid.d_Hyd)
        )
        if self.fluid.MS:
            self.epsilon = (
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

            if self.epsilon > 1e5:
                self.eff_an = 1 - np.exp(-self.tau)
            elif self.epsilon**0.5 < 1e-2 and self.tau < 1 / self.epsilon**0.5:
                self.eff_an = 1 - (1 - self.tau * self.epsilon**0.5) ** 2
            else:
                beta = (1 / self.epsilon + 1) ** 0.5 + np.log(
                    (1 / self.epsilon + 1) ** 0.5 - 1
                )
                max_exp = np.log(np.finfo(np.float64).max)
                beta_tau = beta - self.tau - 1
                if beta_tau > max_exp:
                    print(
                        "Warning: Overflow encountered in exp, input too large.Iterative solver triggered"
                    )
                    # we can use the approximation w=beta_tau-np.log(beta_tau)for the lambert W function but it leads to error up to 40 % in very niche scenarios.

                    def eq(var):
                        cl = var
                        alpha = self.epsilon * self.c_in
                        left = (cl / alpha + 1) ** 0.5 + np.log(
                            (cl / alpha + 1) ** 0.5 - 1 + 1e-10
                        )
                        right = beta - self.tau

                        return abs(left - right)

                    cl = minimize(
                        eq,
                        self.c_in / 2,
                        method="Powell",
                        bounds=[(0, self.c_in)],
                        tol=1e-7,
                    ).x[0]
                    self.eff_an = 1 - (cl / self.c_in)
                else:
                    z = np.exp(beta_tau)
                    w = lambertw(z, tol=1e-10)
                    self.eff_an = 1 - self.epsilon * (w**2 + 2 * w)
                    if self.eff_an.imag != 0:
                        raise ValueError("self.eff_an has a non-zero imaginary part")
                    else:
                        self.eff_an = self.eff_an.real  # get rid of 0*j

            # max_exp = np.log(np.finfo(np.float64).max)
            # beta_tau = beta - self.tau - 1
            # if beta_tau > max_exp:
            #     print("Warning: Overflow encountered in exp, input too large.")
            #     # Handle the overflow case here, e.g., by setting a maximum value
            #     z = np.finfo(np.float64).max
            # else:
            #     z = np.exp(beta_tau)
            # w = lambertw(z, tol=1e-10)
            # self.eff_an = 1 - self.epsilon * (w**2 + 2 * w)
            # if self.eff_an.imag != 0:
            #     raise ValueError("self.eff_an has a non-zero imaginary part")
            # else:
            #     self.eff_an = self.eff_an.real  # get rid of 0*j
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

    def get_flux(self, c: float = None, c_guess: float = 1e-9):
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
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
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
                        * (c / self.fluid.Solubility) ** 0.5
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
                                * (c_wl / self.fluid.Solubility) ** 0.5
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
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
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

                    initial_guess = [(c * 1e-1)]
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
                    self.J_perm = -2 * self.fluid.k_t * c  ## MS factor
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
                            * (self.membrane.K_S * c_ws)
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
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
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
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
                elif self.H < 0.01:
                    # Surface limited V // Mass Transport limited X
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
                            (0, c * (1 + 1e-4)),
                        ],
                        tol=1e-8,
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
                    # Mass transport limited V // Mixed Surface Diffusion  X
                    self.J_perm = -self.fluid.k_t * c  ## LM factor
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

    def get_solid_inventory(self):
        def integrate_c_profile(self):
            r_in = self.fluid.d_Hyd / 2
            r_out = self.fluid.d_Hyd / 2 + self.membrane.thick
            L_min = 0
            L_max = self.geometry.L
            N = 100

            def integrand(r, L):
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
                    c_w = (
                        c * np.exp(L_ch * L) / (dimless2 / self.fluid.k_t + 1)
                    )  # todo check this is liquid conc
                    return (
                        -c_w * np.log(r / r_out) / np.log(r_out / r_in) * 2 * np.pi * r
                    )
                else:
                    tau = 4 * self.fluid.k_t * L / (self.fluid.U0 * self.fluid.d_Hyd)
                    self.epsilon = (
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

                    beta = (1 / self.epsilon + 1) ** 0.5 + np.log(
                        (1 / self.epsilon + 1) ** 0.5 - 1
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
                    conv = (
                        self.c_in / self.fluid.Solubility
                    ) ** 0.5 * self.membrane.K_S
                    c_w_l = alpha * (w**2 + 2 * w) + alpha * (
                        2 - 2 * ((w**2 + 2 * w) + 1) ** 0.5  ## TODO: Check this
                    )
                    return (
                        -np.log(r / r_out)
                        / np.log(r_out / r_in)
                        * 2
                        * np.pi
                        * r
                        * (c_w_l / self.fluid.Solubility) ** 0.5
                        * self.membrane.K_S
                    )

            result, err = integrate.nquad(integrand, [[r_in, r_out], [L_min, L_max]])
            return result

        self.membrane.inv = integrate_c_profile(self)
        return

        # from sympy import Integral, ln, symbols, init_printing, nsolve, exp

        # # Define the symbolic variables
        # r, L = symbols("r L")
        # r_in = self.fluid.d_Hyd / 2
        # r_out = self.fluid.d_Hyd / 2 + self.membrane.thick
        # init_printing(use_unicode=True)
        # if self.fluid.MS == True:
        #     c = (self.c_in / self.fluid.Solubility) ** 0.5 / self.membrane.K_S
        # else:
        #     c = self.c_in / self.fluid.Solubility / self.membrane.K_S

        # # Define the logarithmic function
        # fun1 = -ln(r / r_out) / ln(r / r_in) * c *r * 2 * np.pi

        # # fun2= c*exp(-4*L)
        # fun2 = 1
        # # Perform the integration with respect to r over the interval [r_in, r_out]
        # integral1 = Integral(fun1, (r, r_in, r_out))
        # integral2 = Integral(fun2, (L, 0, self.geometry.L))
        # # Evaluate the integral
        # # integral_value=nsolve(integral1,r, r_out)-nsolve(integral1,r, r_in)

        # print(integral1)

        # integral1.as_sum(5, method="midpoint")
        # integral2.as_sum(5, method="midpoint")
        # result = float(integral1.as_sum(100, method="midpoint")) * float(integral2.as_sum(
        #     100, method="midpoint")
        # )

        # self.membrane.inv = float(result)
        # return

    def get_fluid_inventory(self):
        r_in = self.fluid.d_Hyd / 2

        L_min = 0
        L_max = self.geometry.L
        N = 100

        def integrand(L):
            if self.fluid.k_t is None:
                self.fluid.get_kt(turbulator=self.geometry.turbulator)
            if self.fluid.MS is False:
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

                return self.c_in * np.exp(L_ch * L)
            else:
                tau = 4 * self.fluid.k_t * L / (self.fluid.U0 * self.fluid.d_Hyd)
                self.epsilon = (
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

                beta = (1 / self.epsilon + 1) ** 0.5 + np.log(
                    (1 / self.epsilon + 1) ** 0.5 - 1
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
                                (self.fluid.d_Hyd + 2 * self.membrane.thick)
                                / self.fluid.d_Hyd
                            )
                        )
                    )
                    ** 2
                )
                conv = (self.c_in / self.fluid.Solubility) ** 0.5 * self.membrane.K_S
                c_w_l = alpha * (w**2 + 2 * w)
                return c_w_l

        result, err = integrate.nquad(integrand, [[L_min, L_max]])
        self.fluid.inv = result * np.pi * r_in**2
        return

    def get_inventory(self):
        self.get_solid_inventory()
        self.get_fluid_inventory()
        self.inventory = self.fluid.inv + self.membrane.inv
        return


class Fluid:
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
        self.U0 = U0
        self.k = k
        self.cp = cp
        self.inv = inv

    def update_attribute(
        self, attr_name: str = None, new_value: Union[float, "FluidMaterial"] = None
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

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

    def get_kt(self, turbulator=None):
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
                if turbulator is None:

                    # if Re < 1e4 and Re > 2030:
                    #     Sh = 0.015 * Re**0.83 * Sc**0.42  ## Stempien Thesis pg 155-157 TODO implement different Re ranges
                    if Re > 2030:
                        # Sh = 0.0096 * Re**0.913 * Sc**0.346  ##Getthem paper
                        Sh = 0.023 * Re**0.8 * Sc**0.33
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
                            raise NotImplementedError(
                                "Twisted Tape not implemented yet"
                            )
                        case "Custom":
                            self.k_t = turbulator.k_t_correlation(
                                Re=Re, Sc=Sc, d_hyd=self.d_Hyd, D=self.D
                            )

            else:
                print("k_t is already defined")
        else:
            print("Hydraulic Diameter is not defined")


class Turbulator:
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

    def update_attribute(
        self, attr_name: str = None, new_value: Union[str, dict] = None
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)


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

    def update_attribute(
        self, attr_name: str = None, new_value: Union[str, dict] = None
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)


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

    def update_attribute(
        self, attr_name: str = None, new_value: Union[str, dict] = None
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

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
        k_t = corr.get_k_from_Sh(Sh=Sh, L=d_hyd, D=self.D)
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


class Membrane:
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

    def update_attribute(
        self, attr_name: str = None, new_value: Union[float, "SolidMaterial"] = None
    ):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)
        if self.D_0 is not None and self.E_d is not None:
            if attr_name == "T" and new_value is not None:
                self.D = self.D_0 * np.exp(-self.E_d / (8.617333262145e-5 * self.T))

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

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


# class GLC_Gas:
#     """
#     GLC_Gas class represents a sweep gas in a GLC (Gas-Liquid Contactor) system.

#     Attributes:
#         G_gas (float): The flow rate of the gas.
#         pg_in (float, optional): The Tritium inlet partial pressure of the gas. Default is 0.
#         p_tot (float, optional): The total pressure of the component. Default is 100000 Pa.
#         kla (float, optional): The total (kl*Area) mass transfer coefficient. Default is 0.
#     """

#     def __init__(
#         self,
#         G_gas: float = None,
#         pg_in: float = 0,
#         p_tot: float = 100000,
#         kla: float = 0,
#     ):
#         """
#         Initializes a new instance of the GLC_Gas class.

#         Args:
#             G_gas (float): The flow rate of the gas.
#             pg_in (float, optional): The Tritium inlet partial pressure of the gas. Default is 0.
#             p_tot (float, optional): The total pressure of the component. Default is 100000 Pa.
#             kla (float, optional): The total (kl*Area) mass transfer coefficient. Default is 0.
#         """
#         self.G_gas = G_gas
#         self.pg_in = pg_in
#         self.p_tot = p_tot
#         self.kla = kla

#     def update_attribute(self, attr_name: str, new_value: float):
#         """
#         Updates the value of the specified attribute.

#         Args:
#             attr_name (str): The name of the attribute to update.
#             new_value: The new value for the attribute.
#         """
#         set_attribute(self, attr_name, new_value)

#     def inspect(self, variable_names=None):
#         """
#         Prints the attributes of the component.
#         """
#         print_class_variables(self, variable_names)


# class GLC(Component):
#     """
#     GLC (Gas-Liquid Contact) class represents a gas-liquid contactor component.

#     Args:
#         H (float): Height of the GLC.
#         R (float): Radius of the GLC.
#         L (float): Characteristic Length of the GLC fluid flow.
#         c_in (float): Inlet concentration of the GLC.
#         eff (float, optional): Efficiency of the GLC. Defaults to None.
#         fluid (Fluid, optional): Fluid object representing the liquid phase. Defaults to None.
#         membrane (Membrane, optional): Membrane object representing the membrane used in the GLC. Not super important. Defaults to None.
#         GLC_gas (GLC_Gas, optional): GLC_Gas object representing the gas phase. Defaults to None.

#     Attributes:
#         H (float): Height of the GLC [m].
#         R (float): Radius of the GLC [m].
#         L (float): Length of the GLC [m].
#         GLC_gas (GLC_Gas): GLC_Gas object representing the gas phase.

#     Methods:
#         get_kla_Ring(): Calculates the mass transfer coefficient (kla) for a Raschig Ring matrix.
#     """

#     def __init__(
#         self,
#         H: float = None,
#         R: float = None,
#         L: float = None,
#         c_in: float = None,
#         eff: float = None,
#         fluid: "Fluid" = None,
#         membrane: "Membrane" = None,
#         GLC_gas: "GLC_Gas" = None,
#     ):
#         """
#         Initializes a new instance of the GLC class.

#         Args:
#             H (float): Height of the GLC.
#             R (float): Radius of the GLC.
#             L (float): Characteristic Length of the GLC fluid flow.
#             c_in (float): Inlet concentration of the GLC.
#             eff (float, optional): Efficiency of the GLC. Defaults to None.
#             fluid (Fluid, optional): Fluid object representing the liquid phase. Defaults to None.
#             membrane (Membrane, optional): Membrane object representing the membrane used in the GLC. Not super important. Defaults to None.
#             GLC_gas (GLC_Gas, optional): GLC_Gas object representing the gas phase. Defaults to None.
#         """
#         super().__init__(c_in, eff, fluid, membrane)
#         self.H = H
#         self.R = R
#         self.L = L
#         self.GLC_gas = GLC_gas

#     def get_kla_Ring(self):
#         """
#         Calculates the mass transfer coefficient (kla) for a Raschig ring matrix.

#         The mass transfer coefficient is calculated based on the Reynolds number (Re) and Schmidt number (Sc) of the fluid,
#         as well as the diameter (d) of the ring.

#         Returns:
#             None
#         """
#         d = 2e-3  # Ring diameter
#         Re = corr.Re(rho=self.fluid.rho, u=self.fluid.U0, L=self.L, mu=self.fluid.mu)
#         Sc = corr.Schmidt(D=self.fluid.D, mu=self.fluid.mu, rho=self.fluid.rho)
#         self.GLC_gas.kla = extractor.corr_packed(
#             Re,
#             Sc,
#             d,
#             rho_L=self.fluid.rho,
#             mu_L=self.fluid.mu,
#             L=self.L,
#             D=self.fluid.D,
#         )

#     def inspect(self, variable_names=None):
#         """
#         Prints the attributes of the component.
#         """
#         print_class_variables(self, variable_names)


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

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

    def update_attribute(self, attr_name: str, new_value: float):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)


class SolidMaterial:
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

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

    def update_attribute(self, attr_name: str, new_value: float):
        """
        Updates the value of the specified attribute.

        Args:
            attr_name (str): The name of the attribute to update.
            new_value: The new value for the attribute.
        """
        set_attribute(self, attr_name, new_value)


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
        name: str = None,
    ):
        self.c_in = c_in
        self.Q = Q
        self.TBR = TBR
        self.T_out = T_out
        self.T_in = T_in
        self.fluid = fluid
        self.name = name

    def plot_component(self):
        fig, ax2 = plt.subplots(1, 1, figsize=(10, 5))
        rectangle = plt.Rectangle(
            (0.2, 0.3), 0.55, 0.4, edgecolor="black", facecolor="green", alpha=0.5
        )
        ax2.add_patch(rectangle)
        # Arrow pointing to the left side of the rectangle
        ax2.arrow(
            0.5, 0.7, 0, 0.1, head_width=0.05, head_length=0.1, fc="black", ec="black"
        )
        # Arrow pointing out of the right side of the rectangle
        ax2.arrow(
            0.5, 0.1, 0.0, 0.1, head_width=0.05, head_length=0.1, fc="black", ec="black"
        )
        ax2.set_aspect("equal")
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        if self.name is None:
            ax2.set_title("Component  ")
        else:
            ax2.set_title(self.name)

        # Add text over the arrows
        ax2.text(
            0.7,
            0.8,
            r"$T_o$=" + str(self.T_out) + " K",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.7,
            0.2,
            r"$T_i$=" + str(self.T_in) + " K",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.3,
            0.2,
            f"$c_i$={self.c_in:.4g} $mol/m^3$",
            color="black",
            ha="center",
            va="center",
        )
        ax2.text(
            0.5, 0.6, f"Q={self.Q/1E6:.3g} MW", color="black", ha="center", va="center"
        )
        ax2.text(
            0.5, 0.4, f"TBR={self.TBR:.3g}", color="black", ha="center", va="center"
        )

        ax2.text(
            0.3,
            0.8,
            rf"$c_o$={self.c_out:.4g}$mol/m^3$",
            color="black",
            ha="center",
            va="center",
        )
        ax2.axis("off")
        # Display the plot
        fig.tight_layout()
        return fig

    def inspect(self, variable_names=None):
        """
        Prints the attributes of the component.
        """
        print_class_variables(self, variable_names)

    def update_attribute(self, attr_name: str, new_value: float):
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

    def connect_to_component(self, component2: Union["Component", "BreedingBlanket"]):
        component2.update_attribute("c_in", self.c_out)

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
