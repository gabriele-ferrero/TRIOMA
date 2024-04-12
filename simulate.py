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
        self.initial_conditions = {name: component.tritium_inventory for name, component in component_map.components.items()}
        self.I_startup = component_map.components['Fueling System'].tritium_inventory
        self.y = [list(self.initial_conditions.values())]  # Initialize y with the initial conditions
        self.components = component_map.components
        self.component_map = component_map
        self.interval = self.final_time / 100


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

    def update_timestep(self, dt):
        """
        Update the time step size.

        Args:
        - dt: New time step size.
        """
        self.dt = dt
        self.n_steps = int(self.final_time / dt)


    def adaptive_timestep(self,y_new, y, t, tol=1e-6, max_dt=100, min_dt=1e-6):
        """
        Perform adaptive time stepping.
        """
        p = 1 # Order of the method
        error = np.linalg.norm(y_new - (y + self.dt * self.f(y_new)))
        # Adjust time step size
        dt_new = self.dt * (tol / error)**p
        if abs(t % self.interval) < 10:
            print(f"Time = {t}, Error = {error}, dt = {self.dt}, dt_new = {dt_new}", end='\r')
        # Compute the new time step size based on the definition of adaptive timestep
        dt_new = min(max_dt, max(min_dt, dt_new))
        self.update_timestep(dt_new)

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
        t = 0
        print(self.y[0])
        while t < self.final_time:
            if abs(t % self.interval) < 10:
                print(f"Percentage completed = {abs(t - self.final_time)/self.final_time * 100:.1f}%", end='\r')
            dydt = self.f(self.y[-1])
            y_new = self.y[-1] + self.dt * dydt
            self.time.append(t)
            for i, component in enumerate(self.components.values()):
                component.update_inventory(y_new[i])
            self.component_map.update_flow_rates()
            self.adaptive_timestep(y_new, self.y[-1], t)  # Update the timestep based on the new and old y values
            t += self.dt
            self.y.append(y_new) # append y_new after updating the time step
        self.y = np.array(self.y)
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
