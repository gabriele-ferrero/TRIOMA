class Port:
    def __init__(self, name, incoming_fraction=1.0):
        """
        Initialize a Port object.

        Parameters:
        - name (str): The name of the port.
        - incoming_fraction (float, optional): The fraction of incoming flow to be assigned to this port. Defaults to 1.0.
        """
        self.name = name
        self.flow_rate = 0
        self.incoming_fraction = incoming_fraction
        
    def set_flow_rate(self, flow_rate):
        """
        Set the flow rate of the port.

        Parameters:
        - flow_rate (float): The flow rate to be set.
        """
        self.flow_rate = flow_rate
    