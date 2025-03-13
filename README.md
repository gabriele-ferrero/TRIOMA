# <span style="color:#1E90FF">T</span><span style="color:#32CD32">R</span><span style="color:#FF6347">I</span>OMA (TRItium Object-oriented and Modular Analysis)

[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/gabriele-ferrero/TRIOMA/blob/main/LICENSE) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gabriele-ferrero/TRIOMA/main) [![CI](https://github.com/gabriele-ferrero/TRIOMA/actions/workflows/main.yml/badge.svg)](https://github.com/gabriele-ferrero/TRIOMA/actions) [![codecov](https://codecov.io/gh/gabriele-ferrero/TRIOMA/branch/main/graph/badge.svg)](https://codecov.io/gh/gabriele-ferrero/TRIOMA) ![Stars](https://img.shields.io/github/stars/gabriele-ferrero/TRIOMA.svg?logo=github&label=Stars&logoColor=white) ![Forks](https://img.shields.io/github/forks/gabriele-ferrero/TRIOMA.svg?logo=github&label=Forks&logoColor=white)
TRIOMA is a python-based tool to help engineers to design outer fuel cycles by giving compact, fast and easy pre-built functions to estimate parameters which are crucial for tritium transport, such as the extraction efficiency, losses and inventories.

## Vision

<p style="text-align: justify;">
The objective is to simplify Outer Fuel Cycle (OFC) analysis for Fusion reactors through object-oriented description and efficient pre-built functions, in a easy to understand open-source code to do a quick preliminary estimate of design parameters for an attainable OFC. With TRIOMA it is easy to estimate for an OFC the tritium extraction efficiency of the extractor, the external losses and the inventory.

 The code allows to build OFCs with Packed Tower extractors, heat exchangers, Breeding Blankets and Permeation Against Vacuum Extractors, and accounts the interaction between components, analyzing the OFC as a whole. Both molten salt breeders and liquid metal breeders are implemented, and outlet partial pressure is accounted. The code is verified against the nodal model present in TRIOMA and against FEM models based in COMSOL.

  In the current state, surface effects are implemented only on the outer side of the pipe. This is based on the assumption that inside molten salt fluids the oxyde layer which limits dissociation and recombinarion is not expected to form. In liquid metals Tritium is in the atomic phase, as for the solid metal, so dissociation and recombination are not expected to happen.
</p>

## Getting started

To run copy the repo in a local folder, and run:

```pip install -r requirements.txt```
or pip install trioma directly with
``` pip install TRIOMA ```

## Documentation

Documentation can be viewed at <https://gabriele-ferrero.github.io/TRIOMA/>

## Tutorials

Get started quickly with TRIOMA tutorials at <https://github.com/gabriele-ferrero/TRIOMA/tree/main/Examples>

## Verification

Verification of TRIOMA analytical methods are coming soon

## Contributions

Contributions are welcome from everyone. The project is a work in progress and may intercur significant and structural changes.

## Citations

TRIOMA takes the same approach and often the same equations from the following papers:
Alberghi, Ciro, et al. ["Development of new analytical tools for tritium transport modelling."](<https://www.sciencedirect.com/science/article/pii/S0920379622000837>
)
Fusion Engineering and Design 177 (2022): 113083.
Humrickhouse, Paul W., and Thomas F. Fuerst. [Tritium transport phenomena in molten-salt reactors.](<https://www.osti.gov/biblio/1777267>) No. INL/EXT-20-59927-Rev000. Idaho National Lab.(INL), Idaho Falls, ID (United States), 2020.  

Fuerst, Thomas F., Chase N. Taylor, and Paul W. Humrickhouse. [Tritium Transport Phenomena in Molten-Salt Reactors: Molten Salt Tritium Transport Experiment Design.](<https://www.osti.gov/biblio/1828384>) No. INL/EXT-21-63108-Rev000. Idaho National Lab.(INL), Idaho Falls, ID (United States), 2021.  

Rader, Jordan D., M. Scott Greenwood, and Paul W. Humrickhouse. ["Verification of modelica-based models with analytical solutions for tritium diffusion."](<https://www.tandfonline.com/doi/full/10.1080/00295450.2018.1431505?casa_token=S0I-kCsS6noAAAAA%3A-52Bra2CN56Zg4p9l-l8XXkZXnT0WvPzDI6q-HrQy3NLDPY76wy-UfHlJwZ51VACCmWw7X13Bi-Luc0>) Nuclear Technology 203.1 (2018): 58-65.  

Urgorri, F. R., et al. ["Theoretical evaluation of the tritium extraction from liquid metal flows through a free surface and through a permeable membrane."](<https://iopscience.iop.org/article/10.1088/1741-4326/acbec7/meta>)
