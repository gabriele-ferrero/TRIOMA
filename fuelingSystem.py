from component import Component

class FuelingSystem(Component):
    def __init__(self, name, N_burn, TBE, **kwargs):
        super().__init__(name, residence_time=1, **kwargs)
        self.N_burn = N_burn
        self.TBE = TBE

    # To plasma
    def get_outflow(self):
        return self.N_burn/self.TBE
