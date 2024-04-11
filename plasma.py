from component import Component

class Plasma(Component):
    def __init__(self, name, N_burn, TBE, **kwargs):
        super().__init__(name, residence_time=1, **kwargs)
        self.N_burn = N_burn
        self.TBE = TBE

    # From fueling system
    def get_inflow(self):
        return self.N_burn / self.TBE

    # Exhaust
    def get_outflow(self):
        return (1 - self.TBE)/self.TBE * self.N_burn
    
    def calculate_inventory_derivative(self):
        inflow = self.get_inflow()
        outflow = self.get_outflow()
        dydt = inflow - outflow  + self.tritium_source
        return dydt
