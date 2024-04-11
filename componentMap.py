class ComponentMap:
    """
    A class that represents a component map, which stores information about components and their connections.

    Attributes:
        components (dict): A dictionary that maps component names to component objects.
        connections (dict): A dictionary that maps component names to a dictionary of port names and their connected components and ports.
    """

    def __init__(self):
        self.components = {}
        self.connections = {}

    def add_component(self, component):
        """
        Adds a component to the component map.

        Args:
            component (Component): The component object to be added.
        """
        self.components[component.name] = component

    def connect_ports(self, component1, port1, component2, port2):
        """
        Connects two ports of different components.

        Args:
            component1 (Component): The first component.
            port1 (Port): The port of the first component.
            component2 (Component): The second component.
            port2 (Port): The port of the second component.
        """
        if component1.name not in self.connections:
            self.connections[component1.name] = {}
        if component2.name not in self.connections:
            self.connections[component2.name] = {}

        self.connections[component1.name][port1.name] = (component2.name, port2.name)
        self.connections[component2.name][port2.name] = (component1.name, port1.name)
        if port1 and port2:
                port1.set_flow_rate(component1.get_outflow())
                port2.set_flow_rate(component1.get_outflow() * port2.incoming_fraction)

    def disconnect_ports(self, component1, port1, component2, port2):
        """
        Disconnects two ports of different components.

        Args:
            component1 (Component): The first component.
            port1 (Port): The port of the first component.
            component2 (Component): The second component.
            port2 (Port): The port of the second component.
        """
        if component1.name in self.connections and port1.name in self.connections[component1.name]:
            del self.connections[component1.name][port1.name]
        if component2.name in self.connections and port2.name in self.connections[component2.name]:
            del self.connections[component2.name][port2.name]

    def get_connected_ports(self, component, port):
        """
        Returns the connected component and port for a given component and port.

        Args:
            component (Component): The component.
            port (Port): The port of the component.

        Returns:
            tuple: A tuple containing the connected component and port.
        """
        if component.name in self.connections and port.name in self.connections[component.name]:
            connected_component_name, connected_port_name = self.connections[component.name][port.name]
            connected_component = self.components[connected_component_name]
            connected_port = connected_component.input_ports[connected_port_name]
            return connected_component, connected_port
        else:
            return None, None
        
    def update_flow_rates(self):
        """
        Updates the flow rates of the ports based on the component's outflow and incoming fraction.
        """
        for component_name, ports in self.connections.items():
            for port_name, (connected_component_name, connected_port_name) in ports.items():
                if port_name in self.components[component_name].output_ports:
                    component = self.components[component_name]
                    port = component.output_ports[port_name]
                    connected_component = self.components[connected_component_name]
                    connected_port = connected_component.input_ports[connected_port_name]
                    port.set_flow_rate(component.get_outflow())
                    connected_port.set_flow_rate(component.get_outflow() * connected_port.incoming_fraction)

    def print_connected_map(self):
        """
        Prints the connected map, showing the connections between components and ports.
        """
        for component_name, ports in self.connections.items():
            print(f"Component: {component_name}")
            for port_name, (connected_component_name, connected_port_name) in ports.items():
                print(f"  Port: {port_name} -> Connected Component: {connected_component_name}, Connected Port: {connected_port_name}")
