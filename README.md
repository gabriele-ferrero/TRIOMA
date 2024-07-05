# TRIOMA (TRItium Object-oriented and Modular Analysis)

## Vision

THe objective is to simplify Outer Fuel Cycle (OFC)analysis for Fusion reactors through object-oriented description and efficient pre-built functions, in a easy to understand open-source code to do a quick preliminary estimate of design parameters for an attainable OFC. In the current state, the analysis is limited to pipes with a tritiated liquid medium and zero Tritium partial pressure on the outside, which is the equivalent boundary condition for a constant gas sweep or in presence of water on the outside. In the current state, surface effects are implemented only on the outer side of the pipe. This is based on the assumption that inside molten salt fluids the oxyde layer which limits dissociation and recombinarion is not expected to form. In liquid metals Tritium is in the atomic phase, as for the solid metal, so dissociation and recombination are not expected to happen.
In its current state, the code is capable enough to simulate heat exchangers, Permeation against vacuum and Breeding blanket, together with their interaction.

## Install

To run copy the repo in a local folder, and run:

```pip install -r requirements.txt```

## Contributions

Contributions are welcome from everyone. The project is a work in progress and may intercur significant and structural changes.

## Citations

Alberghi, Ciro, et al. ["Development of new analytical tools for tritium transport modelling."](<https://www.sciencedirect.com/science/article/pii/S0920379622000837>
)
Fusion Engineering and Design 177 (2022): 113083.
Humrickhouse, Paul W., and Thomas F. Fuerst. [Tritium transport phenomena in molten-salt reactors.](<https://www.osti.gov/biblio/1777267>) No. INL/EXT-20-59927-Rev000. Idaho National Lab.(INL), Idaho Falls, ID (United States), 2020.  

Fuerst, Thomas F., Chase N. Taylor, and Paul W. Humrickhouse. [Tritium Transport Phenomena in Molten-Salt Reactors: Molten Salt Tritium Transport Experiment Design.](<https://www.osti.gov/biblio/1828384>) No. INL/EXT-21-63108-Rev000. Idaho National Lab.(INL), Idaho Falls, ID (United States), 2021.  

Rader, Jordan D., M. Scott Greenwood, and Paul W. Humrickhouse. ["Verification of modelica-based models with analytical solutions for tritium diffusion."](<https://www.tandfonline.com/doi/full/10.1080/00295450.2018.1431505?casa_token=S0I-kCsS6noAAAAA%3A-52Bra2CN56Zg4p9l-l8XXkZXnT0WvPzDI6q-HrQy3NLDPY76wy-UfHlJwZ51VACCmWw7X13Bi-Luc0>) Nuclear Technology 203.1 (2018): 58-65.  

Urgorri, F. R., et al. ["Theoretical evaluation of the tritium extraction from liquid metal flows through a free surface and through a permeable membrane."](<https://iopscience.iop.org/article/10.1088/1741-4326/acbec7/meta>)
