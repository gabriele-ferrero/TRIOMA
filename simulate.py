import numpy as np

class Simulate:
    def __init__(self, dt, final_time, component_map):
        """
        Initialize the Simulate class.

        Args:
        - dt: Time step size.
        - final_time: Final simulation time.
        - component_map: Mapping of component names to Component objects.
        """
        self.dt = dt
        self.final_time = final_time
        self.time = []
        self.n_steps = int(final_time / dt)
        self.initial_conditions = {name: component.tritium_inventory for name, component in component_map.components.items()}
        self.I_startup = component_map.components['Fueling System'].tritium_inventory
        self.y = np.zeros((self.n_steps + 1, len(self.initial_conditions)))
        self.components = component_map.components
        self.component_map = component_map

    def run(self):
        """
        Run the simulation.

        Returns:
        - t: Array of time values.
        - y: Array of component inventory values.
        """
        while True:
            self.y[0] = [component.tritium_inventory for component in self.components.values()] # self.initial_conditions, possibly updated by the restart method
            t,y = self.forward_euler()
            if self.components['Fueling System'].tritium_inventory < 0:
                print("Error: Tritium inventory in Fueling System is below zero.")
                self.update_I_startup()
                self.restart()
            else:
                return t,y
            
    def update_I_startup(self):
        """
        Update the initial tritium inventory of the Fueling System component.
        """
        self.I_startup += 0.1
        self.initial_conditions['Fueling System'] = self.I_startup

    def restart(self):
        """
        Restart the simulation by resetting time and component inventory.
        """
        self.time = []
        self.y = np.zeros((self.n_steps + 1, len(self.initial_conditions)))
        for component, initial_condition in zip(self.components.values(), self.initial_conditions.values()):
            component.tritium_inventory = initial_condition

    def forward_euler(self):
        """
        Perform the forward Euler integration method.

        Returns:
        - time: Array of time values.
        - y: Array of component inventory values.
        """
        print(self.y[0])
        for n in range(self.n_steps):
            t = n * self.dt
            if n % 1000 == 0:
                print(f"Percentage completed = {n/self.n_steps * 100:.1f}%", end='\r')
            dydt = self.f(self.y[n])
            self.y[n+1] = self.y[n] + self.dt * dydt
            self.time.append(t)
            for i, component in enumerate(self.components.values()):
                component.update_inventory(self.y[n+1][i])
            self.component_map.update_flow_rates()
        return [self.time, self.y[:-1,:]]


    def f(self, y):
        """
        Calculate the derivative of component inventory.

        Args:
        - y: Array of component inventory values.

        Returns:
        - dydt: Array of derivative values.
        """
        dydt = np.zeros_like(y)
        for i, component in enumerate(self.components.values()):
            dydt[i] += component.calculate_inventory_derivative()
        return dydt
