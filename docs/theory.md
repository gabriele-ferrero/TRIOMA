# Theory

## Blanket

The Breeding blanket is the most important component for an OFC, and is the source term for tritium in the loop. However, TRIOMA can also simulate experimental layouts, therefore it is not a mandatory component.

The breeding blanket is also the most complex component for the OFC, and TRIOMA can't handle all the complexities derived by geometries and pulsed operations. This is mostly why steady state is considered, because the BB can be modeled as a balance of mass, hydrogen, and heat.
The necessary parameters to describe the breeding blanket are:

- $c_{in}$
- Tritium Breeding ratio
- Plasma power
- mass flow rate
- Fluid density

from plasma power, the tritons per second of the plasma are evaluated as neutron per seconds $\cdot TBR$:
$t_{gen}=\frac{Q_{plasma}}{17.6\left[MeV\right] \cdot1,6022\cdot 10^{-19}\left[\frac{J}{eV}\right]}\cdot \frac{TBR}{N_A}  $
in $ mol_T/s$

Then $c_{out}$ is evaluated as:

$c_{out,lm}[mol_{T}/m^3]=c_{in}+\frac{t_{gen}}{\dot{m}_{breeder}/\rho_{breeder}}$

for liquid metals, where hydrogen isotopes are in the atomic form T. In molten salts, the molecular form $T_2$ is found, therefore is accounted in the balance of tritons per second:

$ c_{out,ms} [mol_{T_2}/m^3]=c_{in}+\frac{t_{gen}}{2\cdot \dot{m}_{breeder}/\rho_{breeder}}$

For the Liquid Immersion Blanket of ARC-class reactors, the tritium carrier is also the power exhaust fluid. Therefore, it is possible to evaluate the mass flow rate from the BreedingBlanket.get_flowrate() method with a given $c_p$ and $\Delta T$.

An application of the BreedingBlanket class can be seen at the circuit class [tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/circuit_tutorial.ipynb)

## Component

Component is the name used in TRIOMA for pipe components. This Class can describe the behaviour of heat exchangers, connection pipes and Permeation Against Vacuum Extractors, making it very useful for OFCs modeling.
In the pipe fluid, hydrogen transport is regulated by convective transport and is described as  

$
J = k_t \cdot (c_{bulk} - c_{w,l})
$

where $k_t$ is the mass transfer coefficient in $m/s$, $c_{bulk}$ is the hydrogen concentration in the center of the liquid and $c_{w,l}$ is the concentration at the wall.
Then, diffusive effects in the pipe wall are described through Richardson Law:  

 $J=\frac{\phi}{r \cdot \ln\left(\frac{r_o}{r}\right)} \cdot \left( {\sqrt{p_{in,Q_2}}} - \sqrt{p_{out,Q_2}} \right)$

 Where $\phi$ is the material permeability ($\phi=D\cdot K_S $), and ${p_{in,Q_2}}$ is the hydrogen partial pressure inside the pipe wall interface. In steady state condition, the two fluxes are equal. However, ${p_{in,Q_2}}$ is connected with $c_{w,l}$ with two different functions depending by the fluid.
 If the fluid is a liquid metal:

 $\sqrt{p_{Q_2}}=\frac{c_{w,l}}{K_{S,l}}=\frac{c_{w,s}}{K_{S,s}}$

 being $K_S$ the Sieverts solubility constant, $K_{S,l}$ in the liquid and $K_{S,s}$ in the solid.
 For a molten salt:

 $    c_{Q_2}={p_{Q_2}}{K_H}$

 and

$ p_{Q_2}=\frac{c_{Q_2,l}}{K_H}=\left(\frac{c_{Q,s}}{K_S}\right)^2$

 being $K_H$ Henry Solubility constant.
 In a fluid with hydrogen content, the change in concentration is described as:

 $
  \frac{\delta c_{bulk}(z)}{\delta z}=\frac{\pi d J}{v \pi \frac{d^2}{4}}=   \frac{-4J}{vd}
  $

Then the solution for J is substituted in the previous equation, to evaluate $c_{out}$ and $\eta$.
To do so, the component class needs all the ingriedients which correctly describe the component, which in TRIOMA are the Fluid, the Membrane and the Geometry classes. Without all the inputs, the hydrogen transport obviously cannot be solved.

### Membrane class

The membrane class represents the pipe wall properties. The pipe wall characteristics mandatory for hydrogen transport are:

- Thickness [t]
- Diffusivity [D]
- Solubility [K_S]
- Operating temperature [T]

### Fluid class

The fluid class represents the liquid carrying the hydrogen isotope. The fluid characteristics mandatory for hydrogen transport are:

- Operating temperature [T]
- Solubility in the liquid [Solubility]
- Diffusivity in the liquid [D]
- Hydraulic diameter [d_Hyd]
- Density [rho]
- Velocity [U_0]
- Viscosity [mu]
- MS: is the fluid a molten salt (True) of a liquid metal (False)

By these properties, $k_t$ is evaluated through correlations

### Geometry class

The geometry class represents the geometry of the component. It repeats some variables from fluid and membrane classes, as well as some of components.
The geometry characteristics needed for hydrogen transport are:

- Length of the pipe [L]
- Inner diameter [D]
- Thickness [t]
- number of pipes [n_pipes]

### Component class resume

The component class represents the pipe.
The component characteristics needed for hydrogen transport are:

- Geometry
- Fluid
- Membrane
- $c_{in}$ [c_in]

An application of TRIOMA components is shown at the [tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/tutorial.ipynb), the recirculation [tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/recirculation.ipynb)

### Note on Heat Exchangers

Heat exchangers have a temperature profile, which is not represented in TRIOMA 0D objects. Therefore, first the component is defined, then the split_HX function is used to build a Circuit made of shorter pipes, each with a temperature evaluated through an HX dimensioning built in TRIOMA.
To know how many segments are needed to evaluate correctly the transport in the HX, run the Component.converge_slit_HX() method. Example is seen in the [HX Tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/Heat_exchanger_tutorial.ipynb)

![Temperature distribution in HX](Figures\HX_temperature_profile.png)

*Figure 1: Temperature distrubution in a HX component.*

## Circuit

The circuit class is used to connect components and evaluate the global behaviour of the OFC considering interactions between components.
To simulate a circuit it is necessary to include

- components : a list of components( in order) which can be Component classes, GLC classes, Circuit classes (it extracts all its components and appends them to the list), Breeding blanket classes.
- closed: a boolean variable to decide if the last component of the list must be connected to the first and make a closed loop.

In particular, only one Breeding blanket is accepted. One real BB is already hard enough to build.

The connection simply is:
$c_{in}^{n+1}=c_{out}^n$
for the n component and the following (n+1) component.

The closed boolean makes the circuit run until the steady state solution is found , in which the hydrogen introduced in the circuit and the hydrogen extracted are the same quantity. In particular it is evaluated as :

$c_{in,BB}=c_{out,OFC}(c_{in,BB})$

where $c_{out,OFC}(c_{in,BB})$ is the outlet concentration of the OFC when (c_{in,BB}) is the inlet at the breeding blanket.

If the circuit is open (closed=False), the inlet concentration of the first component overrides all other inlet concentration through the connection function.
An application of the Circuit class is shown at the circuit [tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/circuit_tutorial.ipynb).

## GLC

Gas liquid contactors are packed tower extractors in which a gas is put directly in contact with the fluid, and a packing increases mass transfer coefficients between the two. The modeling method is the Height of transfer Units and Number of transfer units (HTU-NTU) method.

![Sketch for a GLC](Figures\GLC.png)
*Figure 2:Sketch for a GLC component.*

Mass transfer in GLC is described by the mass transfer coefficient $k_l$ multiplied by contact area a in  $m^2/m^3$
In liquid metals, due to different forms , the Q flux is evaluated as:
$  0.5\cdot u_gdc_g=-u_Ldc_L=k_La(c_L-c_L^*)dz $
$c_L^*$ is the liquid concentration when in equilibrium with the partial pressure of the sweep gas.

$\left(\frac{c_L^*}{K_S}\right)^2=p_eq=c_GRT$

Then combining and integrating the two equations gives:

$ NTU_{LM}=\int^{c_{L,in}}_{c_{L,out}} \frac{dc_L}{c_L-K_S\cdot((\frac{u_L}{2u_G}RT)^{1/2}-(c_L-c_{L,out}+c_{g,in}\frac{u_G}{2\cdot u_L})^{1/2}}$
 and :

$ HTU=\frac{u_L}{k_La}$

 $Z=HTU\cdot NTU$

 Where Z is the height of the packing. Therefore TRIOMA combines these equations in different ways:

- If Z is known, and $c_{out}$ are known, $k_La$ can be evaluated through GLC.get_kla_from_cout().
- If $k_La$ and Z are known, $c_{out}$ is evaluated with get_c_out().
- If $c_{out}$ and $k_La$ are known, the method get_z_from_eff() can evaluate the needed height to get a certain efficiency therefore the desired $c_{out}$.

For molten salt systems, the functions are the same but the definition of $c_L^*$ therefore NTU is different.

$NTU_{MS}=\int^{c_{L,in}}_{c_{L,out}} \frac{dc_L}{c_L-K_H\frac{u_L}{u_G}RT\cdot\left(c_L-c_{L,out}+c_{g,in}\frac{u_G}{u_L}\right)}$

To simulate a GLC, the following variables must be defined:

- GLC_GAS
- Fluid (for Solubility, MS boolean)
- two variables out of Z, k_La

An application for GLC extractor is shown at the GLC [tutorial](https://github.com/gabriele-ferrero/TRIOMA/blob/main/Examples/GLC_tutorial.ipynb)

### GLC_GAS class

This class is needed in the GLC class to store characteristics of the sweep gas. To simulate a GLC, the following variables are needed:

- Gas flowrate [G_gas]
- inlet hydrogen isotope partial pressure [pg_in]
- gas pressure [p_tot]

from the hydrogen balance, pg_out is evaluated when using the get_c_out() of GLC class.
