from TRIOMA.tools.TriomaClass import TriomaClass
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
from TRIOMA.tools.Extractors.GasLiquidContactor import GLC_Gas, GLC

import TRIOMA.tools.correlations as corr
import matplotlib.pyplot as plt

from typing import Union
import matplotlib.pyplot as plt
from scipy import integrate
from TRIOMA.tools.Extractors import extractor
from TRIOMA.tools.Extractors.PAV import Component
from TRIOMA.tools.BreedingBlanket import BreedingBlanket
from TRIOMA.tools.Circuit import Circuit
import numpy as np


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


def connect_to_component(
    self, component2: Union["Component", "BreedingBlanket", "GLC"] = None
):
    """sets the inlet conc of the object component equal to the outlet of self"""
    component2.update_attribute("c_in", self.c_out)


Component.converge_split_HX = converge_split_HX
Component.split_HX = split_HX
Component.connect_to_component = connect_to_component
BreedingBlanket.connect_to_component = connect_to_component
GLC.connect_to_component = connect_to_component
