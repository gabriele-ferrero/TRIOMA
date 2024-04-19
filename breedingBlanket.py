from component import Component

class BreedingBlanket(Component):
    def __init__(self, name, residence_time,  N_burn, TBR, initial_inventory=0, non_radioactive_loss=0.0001, **kwargs):
        super().__init__(name, residence_time, initial_inventory, non_radioactive_loss, **kwargs)
        self.N_burn = N_burn
        self._TBR = TBR  # Initialize _TBR directly
        self.tritium_source = self.N_burn * self.TBR

    @property
    def TBR(self):
        return self._TBR

    @TBR.setter
    def TBR(self, value):
        self._TBR = value
        self.tritium_source = self.N_burn * self.TBR


