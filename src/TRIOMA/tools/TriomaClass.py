from typing import Optional


class TriomaClass:
    """Base class for TRIOMA components with inspection and attribute management."""

    def inspect(self, variable_names: Optional[str] = None, tab: int = 0) -> None:
        """
        Print specified variables of the class.

        If a variable is a class itself, calls this function recursively.
        If variable_names is None, prints all variables.

        Args:
            variable_names: Name of the variable to print. Prints all if None.
            tab: Indentation level for nested inspection.
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

    def update_attribute(self, attr_name: str, new_value: float) -> None:
        """
        Set the specified attribute to a new value.

        Searches recursively through this object and all nested TriomaClass objects.

        Args:
            attr_name: The name of the attribute to set.
            new_value: The new value for the attribute.

        Raises:
            ValueError: If the attribute doesn't exist in this object or nested objects.
        """
        if not self._update_attribute_recursive(attr_name, new_value):
            raise ValueError(f"'{attr_name}' is not an attribute of {self.__class__.__name__}")

    def _update_attribute_recursive(self, attr_name: str, new_value: float) -> bool:
        """
        Recursively search for and update an attribute in this object or nested objects.

        Args:
            attr_name: The attribute name to search for.
            new_value: The new value to set.

        Returns:
            True if attribute was found and updated, False otherwise.
        """
        # Check current object
        if hasattr(self, attr_name):
            setattr(self, attr_name, new_value)

            # For special attributes, also update in nested objects
            if attr_name == "n_pipes":
                for nested_obj in self.__dict__.values():
                    if isinstance(nested_obj, TriomaClass):
                        nested_obj._update_attribute_recursive(attr_name, new_value)

            return True

        # Search nested TriomaClass objects
        for nested_obj in self.__dict__.values():
            if isinstance(nested_obj, TriomaClass):
                if nested_obj._update_attribute_recursive(attr_name, new_value):
                    return True

        return False
