from component import Component

class FuelingSystem(Component):
    def __init__(self, name, N_burn, TBE, **kwargs):
        """
        Initialize a FuelingSystem object.

        Args:
            name (str): The name of the fueling system.
            N_burn (float): The number of burn cycles.
            TBE (float): The time between each burn cycle.
            **kwargs: Additional keyword arguments.

        Returns:
            None
        """
        super().__init__(name, residence_time=1, **kwargs)
        self.N_burn = N_burn
        self.TBE = TBE

    # To plasma
    def get_outflow(self):
        """
        Calculate the outflow rate of the fueling system.

        Returns:
            float: The outflow rate.
        """
        return self.N_burn/self.TBE
