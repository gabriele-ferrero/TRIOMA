
from fuelingSystem import FuelingSystem
from component import Component
from plasma import Plasma
from componentMap import ComponentMap
from matplotlib import pyplot as plt
from simulate import Simulate
from tools.utils import visualize_connections

LAMBDA = 1.73e-9 # Decay constant for tritium
N_burn = 9.3e-7 # Tritium burn rate in the plasma
TBR = 1.1
tau_ofc = 2 * 3600
tau_ifc = 5 * 3600
tau_tes = 24 * 3600
tau_HX = 1 * 3600
I_startup = 0.7 
TBE = 0.02
tes_efficiency = 0.9

# Define components
fueling_system = FuelingSystem("Fueling System", N_burn, TBE, initial_inventory=I_startup)
BB = Component("BB", tau_ofc, initial_inventory=0, tritium_source=N_burn * TBR)
IFC = Component("IFC", tau_ifc)
plasma = Plasma("Plasma", N_burn, TBE) 
TES = Component("TES", residence_time = tau_tes)
HX = Component("HX", residence_time = tau_HX)

# Define ports
port1 = fueling_system.add_output_port("Fueling to Plasma")
port2 = plasma.add_input_port("Port 2")
port3 = plasma.add_output_port("Plasma to IFC")
port4 = IFC.add_input_port("Port 4")
port5 = IFC.add_output_port("IFC to Fueling System")
port6 = BB.add_output_port("OFC to TES")
port7 = fueling_system.add_input_port("Port 7")
port8 = IFC.add_input_port("Port 8")
port9 = TES.add_output_port("TES to Fueling System")
port10 = TES.add_output_port("TES to HX")
port11 = TES.add_input_port("Port 11")
port12 = fueling_system.add_input_port("Port 12", incoming_fraction=tes_efficiency)
port13 = HX.add_input_port("Port 13", incoming_fraction=1-tes_efficiency)
port14 = HX.add_output_port("Port 14")
port15 = BB.add_input_port("Port 15")

# Add components to component map
component_map = ComponentMap()
component_map.add_component(fueling_system)
component_map.add_component(BB)
component_map.add_component(IFC)
component_map.add_component(plasma)
component_map.add_component(TES)
component_map.add_component(HX)

# Connect ports
component_map.connect_ports(fueling_system, port1, plasma, port2)
component_map.connect_ports(plasma, port3, IFC, port4)
component_map.connect_ports(IFC, port5, fueling_system, port7)
component_map.connect_ports(BB, port6, TES, port11)
component_map.connect_ports(TES, port9, fueling_system, port12)
component_map.connect_ports(TES, port10, HX, port13)
component_map.connect_ports(HX, port14, BB, port15)

component_map.print_connected_map()
visualize_connections(component_map)

simulation = Simulate(0.1, 1e5, component_map)
t, y = simulation.run()
plt.figure()
plt.plot(t, y)
plt.legend(component_map.components.keys())
print(f"Component inventories {y[-1]}")