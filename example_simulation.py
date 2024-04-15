from fuelingSystem import FuelingSystem
from component import Component
from plasma import Plasma
from componentMap import ComponentMap
from matplotlib import pyplot as plt
from simulate import Simulate
from utils import visualize_connections

LAMBDA = 1.73e-9 # Decay constant for tritium
N_burn = 9.3e-7 # Tritium burn rate in the plasma
TBR = 1.1
tau_ofc = 2 * 3600
tau_ifc = 4 * 3600
I_startup = 0.7 
TBE = 0.02


component1 = FuelingSystem("Fueling System", N_burn, TBE, initial_inventory=I_startup)
component2 = Component("OFC", tau_ofc, initial_inventory=0, tritium_source=N_burn * TBR)
component3 = Component("IFC", tau_ifc)
component4 = Plasma("Plasma", N_burn, TBE)
port1 = component1.add_output_port("Port 1")
port2 = component4.add_input_port("Port 2")
port3 = component4.add_output_port("Port 3")
port4 = component3.add_input_port("Port 4")
port5 = component3.add_output_port("Port 5")
port6 = component2.add_output_port("Port 6")
port7 = component1.add_input_port("Port 7")
port8 = component3.add_input_port("Port 8")

component_map = ComponentMap()
component_map.add_component(component1)
component_map.add_component(component2)
component_map.add_component(component3)
component_map.add_component(component4)
component_map.connect_ports(component1, port1, component4, port2)
component_map.connect_ports(component4, port3, component3, port4)
component_map.connect_ports(component3, port5, component1, port7)
component_map.connect_ports(component2, port6, component3, port8)
component_map.print_connected_map()
visualize_connections(component_map)

simulation = Simulate(0.1, 1e5, component_map)
t, y = simulation.run()
# plt.figure()
# plt.plot(t, y)
# plt.legend(component_map.components.keys())
# print(y[-1,:])