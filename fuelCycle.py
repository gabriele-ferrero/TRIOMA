
from fuelingSystem import FuelingSystem
from component import Component
from plasma import Plasma
from breedingBlanket import BreedingBlanket
from componentMap import ComponentMap
from matplotlib import pyplot as plt
from simulate import Simulate
from tools.utils import visualize_connections

LAMBDA = 1.73e-9 # Decay constant for tritium
AF = 0.7
N_burn = 9.3e-7 * AF # Tritium burn rate in the plasma
TBR = 1.01
tau_ofc = 2 * 3600
tau_ifc = 5 * 3600
tau_tes = 24 * 3600
tau_HX = 1 * 3600
tau_FW = 1000
tau_div = 1000
tau_ds = 3600
I_startup = 1.5
TBE = 0.02
tes_efficiency = 0.9
final_time = 2.1 * 3600 * 24 * 365 # NB: longer than doubling time
hx_to_fw = 0.33
hx_to_div = 0.33
hx_to_ds = 1e-4
hx_to_BB = 1 - hx_to_fw - hx_to_div - hx_to_ds

q = 0.25
t_res = 24 * 3600
I_reserve = N_burn / TBE * q * t_res


# Define components
fueling_system = FuelingSystem("Fueling System", N_burn, TBE, initial_inventory=I_startup)
BB = BreedingBlanket("BB", tau_ofc, initial_inventory=0, N_burn = N_burn, TBR = TBR)
FW = Component("FW", residence_time = 600)
divertor = Component("Divertor", residence_time = 600)
IFC = Component("IFC", tau_ifc)
plasma = Plasma("Plasma", N_burn, TBE) 
TES = Component("TES", residence_time = tau_tes)
HX = Component("HX", residence_time = tau_HX)
DS = Component("DS", residence_time = tau_ds)

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
port14 = HX.add_output_port("HX to BB")
port15 = BB.add_input_port("Port 15", incoming_fraction= hx_to_BB)
port16 = FW.add_input_port("Port 16", incoming_fraction=hx_to_fw)
port17 = FW.add_output_port("FW to BB")
port18 = divertor.add_input_port("Port 18", incoming_fraction=hx_to_div)
port19 = divertor.add_output_port("Divertor to FW")
port20 = HX.add_output_port("HX to FW")
port21 = HX.add_output_port("HX to div")
port22 = BB.add_input_port("Port 22")
port23 = BB.add_input_port("Port 23")
port24 = DS.add_input_port("Port 24", incoming_fraction=hx_to_ds)
port25 = DS.add_output_port("DS to IFC")
port26 = HX.add_output_port("HX to DS")
port27 = IFC.add_input_port("Port 27")  

# Add components to component map
component_map = ComponentMap()
component_map.add_component(fueling_system)
component_map.add_component(BB)
component_map.add_component(IFC)
component_map.add_component(plasma)
component_map.add_component(TES)
component_map.add_component(HX)
component_map.add_component(FW)
component_map.add_component(divertor)
component_map.add_component(DS)

# Connect ports
component_map.connect_ports(fueling_system, port1, plasma, port2)
component_map.connect_ports(plasma, port3, IFC, port4)
component_map.connect_ports(IFC, port5, fueling_system, port7)
component_map.connect_ports(BB, port6, TES, port11)
component_map.connect_ports(TES, port9, fueling_system, port12)
component_map.connect_ports(TES, port10, HX, port13)
component_map.connect_ports(HX, port14, BB, port15)
component_map.connect_ports(HX, port20, FW, port16)
component_map.connect_ports(HX, port21, divertor, port18)
component_map.connect_ports(FW, port17, BB, port22)
component_map.connect_ports(divertor, port19, BB, port23)
component_map.connect_ports(HX, port26, DS, port24)
component_map.connect_ports(DS, port25, IFC, port27)

component_map.print_connected_map()
visualize_connections(component_map)
print(f'Startup inventory is: {fueling_system.tritium_inventory}')
simulation = Simulate(dt=0.1, final_time=final_time, I_reserve=I_reserve, component_map=component_map)
t, y = simulation.run()
plt.figure()
plt.plot(t, y)
plt.legend(component_map.components.keys())
print(f"Component inventories {y[-1]}")