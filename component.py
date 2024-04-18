from port import Port
LAMBDA = 1.73e-9 # Decay constant for tritium

class Component:
    """
    Represents a component in a fuel cycle system.
    """

    def __init__(self, name, residence_time, initial_inventory=0, tritium_source=0, non_radioactive_loss=1e-4):
        """
        Initializes a Component object.

        Args:
            name (str): The name of the component.
            residence_time (float): The residence time of the component in seconds.
            initial_inventory (float, optional): The initial tritium inventory of the component. Defaults to 0.
            tritium_source (float, optional): The tritium source rate for the component. Defaults to 0.
        """
        self.name = name
        self.residence_time = residence_time
        self.input_ports = {}  # Dictionary where the key is the port name and the value is the port object
        self.output_ports = {}  # Dictionary where the key is the port name and the value is the port object
        self.tritium_inventory = initial_inventory
        self.tritium_source = tritium_source
        self.non_radioactive_loss = non_radioactive_loss

    def add_input_port(self, port_name, incoming_fraction=1.0):
        """
        Adds an input port to the component.

        Args:
            port_name (str): The name of the input port.
            incoming_fraction (float, optional): The fraction of incoming flow to the port. Defaults to 1.0.

        Returns:
            Port: The created input port object.
        """
        if not (0 <= incoming_fraction <= 1):
            raise ValueError("Incoming fraction must be between 0 and 1")
        port = Port(port_name)
        port.incoming_fraction = incoming_fraction
        self.input_ports[port_name] = port
        return port

    def add_output_port(self, port_name):
        """
        Adds an output port to the component.

        Args:
            port_name (str): The name of the output port.

        Returns:
            Port: The created output port object.
        """
        port = Port(port_name)
        self.output_ports[port_name] = port
        return port

    def __str__(self):
        """
        Returns a string representation of the component.

        Returns:
            str: The string representation of the component.
        """
        return f"{self.name}: Residence Time = {self.residence_time}, Tritium Inventory = {self.tritium_inventory}"

    def add_tritium(self, amount):
        """
        Adds tritium to the component's inventory.

        Args:
            amount (float): The amount of tritium to add.
        """
        self.tritium_inventory += amount

    def remove_tritium(self, amount):
        """
        Removes tritium from the component's inventory.

        Args:
            amount (float): The amount of tritium to remove.

        Returns:
            float: The actual amount of tritium removed.
        """
        if self.tritium_inventory >= amount:
            self.tritium_inventory -= amount
            return amount
        else:
            removed_amount = self.tritium_inventory
            self.tritium_inventory = 0
            return removed_amount
        
    def get_inflow(self):
        """
        Calculates the total inflow rate to the component.

        Returns:
            float: The total inflow rate.
        """
        inflow = 0
        for port in self.input_ports.values():
            inflow += port.flow_rate
        return inflow   
        
    def get_outflow(self):
        """
        Calculates the outflow rate from the component.

        Returns:
            float: The outflow rate.
        """
        return self.tritium_inventory / self.residence_time
    
    def calculate_inventory_derivative(self):
        """
        Calculates the derivative of the tritium inventory with respect to time.

        Returns:
            float: The derivative of the tritium inventory.
        """
        inflow = self.get_inflow()
        outflow = self.get_outflow()
        decay = self.tritium_inventory * LAMBDA 
        dydt = inflow - outflow * (1 + self.non_radioactive_loss) - decay + self.tritium_source
        return dydt
    
    def update_inventory(self, new_value):
        """
        Updates the tritium inventory of the component.

        Args:
            new_value (float): The new value of the tritium inventory.
        """
        self.tritium_inventory = new_value