from TRIOMA.tools.Extractors.PAV import Component
from TRIOMA.tools.BreedingBlanket import BreedingBlanket
from TRIOMA.tools.Extractors.GasLiquidContactor import GLC
from TRIOMA.tools.TriomaClass import TriomaClass
import matplotlib.pyplot as plt
from typing import Union


class Circuit(TriomaClass):
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
                if isinstance(element, Union[Component, BreedingBlanket, GLC]):
                    vec_components.append(element)
                elif isinstance(element, Circuit):
                    for comp in element.components:
                        vec_components.append(comp)
                else:
                    raise ValueError("Invalid component type")
        self.components = vec_components
        self.closed = closed

    def add_component(
        self, component: Union["Component", "BreedingBlanket", "Circuit", "GLC"]
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

    def get_circuit_pumping_power(self):
        """
        Calculates the pumping power required for the circuit.

        Returns:
            float: The pumping power required for the circuit in W.
        """
        pumping_power = 0
        for component in self.components:
            if isinstance(component, Component):
                component.get_pressure_drop()
                component.get_pumping_power()
                pumping_power += component.pumping_power
        self.pumping_power = pumping_power
        return self.pumping_power

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
            if isinstance(component, GLC):
                component.get_c_out()
                if i != len(self.components) - 1:
                    component.connect_to_component(self.components[i + 1])
            if isinstance(component, Component):
                component.use_analytical_efficiency(p_out=component.p_out)
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

    def get_inventory(self, flag_an=True):
        """
        Calculates the inventory (in mol) of the circuit based on the components present.

        Returns:
            circuit.inv (float): The inventory of the circuit.

        Raises:
            None

        Example Usage:
            circuit.get_circuit_inventory()

        """
        inventory = 0
        for component in self.components:
            if isinstance(component, Component):
                component.get_inventory(flag_an=flag_an)
                inventory += component.inv
        self.inv = inventory

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
                if isinstance(component, GLC):
                    component.get_c_out()
                    if i != len(self.components) - 1:
                        component.connect_to_component(self.components[i + 1])
                if isinstance(component, Component):
                    component.use_analytical_efficiency(p_out=component.p_out)
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

    def estimate_cost(self, metal_cost=0, fluid_cost=0):
        """
        Estimates the cost of the circuit.
        metal_cost: cost of the metal in $/m^3 in a vector ordered as the components
        fluid_cost: cost of the fluid in $/m^3 in a vector ordered as the components
        returns the cost of the component
        """
        cost_circuit = 0
        for i, component in enumerate(self.components):
            if isinstance(component, Component):
                component.estimate_cost(
                    metal_cost=metal_cost[i], fluid_cost=fluid_cost[i]
                )
                cost_circuit += component.cost
        self.cost = cost_circuit
        return self.cost
