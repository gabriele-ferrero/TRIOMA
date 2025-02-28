class TriomaClass:

    def inspect(self, variable_names=None, tab: int = 0):
        """
        Prints specified variables of a class. If a variable is a class itself,
        calls this function recursively for the internal class. If variable_names is None,
        prints all variables.

        Args:
            instance: The class instance.
            variable_names (list of str, optional): Names of the variables to print. Prints all if None.
        """
        indent = "    " * tab
        for attr_name, attr_value in self.__dict__.items():
            if variable_names is None or attr_name.lower() == variable_names.lower():
                if type(attr_value) is TriomaClass:
                    tab += 1
                    print(
                        f"{indent}{attr_name} is a {type(attr_value)} class, printing its variables:"
                    )
                    TriomaClass.inspect(attr_value, variable_names, tab=tab)
                    tab -= 1
                else:
                    print(f"{indent}{attr_name}: {attr_value}")

    def update_attribute(self, attr_name: str, new_value: float):
        """
        Sets the specified attribute to a new value.

        Args:
            attr_name (str): The name of the attribute to set.
            new_value: The new value for the attribute.
        """
        if hasattr(self, attr_name):
            setattr(self, attr_name, new_value)
            if attr_name == "n_pipes":
                for attr, value in self.__dict__.items():
                    if isinstance(value, object) and hasattr(value, attr_name):

                        setattr(value, attr_name, new_value)
            return
        else:
            for attr, value in self.__dict__.items():
                if isinstance(value, object) and hasattr(value, attr_name):
                    setattr(value, attr_name, new_value)
                    return
        raise ValueError(
            f"'{attr_name}' is not an attribute of {self.__class__.__name__}"
        )
