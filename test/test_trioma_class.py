import unittest
from io import StringIO
from unittest.mock import patch

from TRIOMA.tools.TriomaClass import TriomaClass


class TestInspect(unittest.TestCase):
    """Essential tests for TriomaClass.inspect() method."""

    def test_inspect_all_attributes(self):
        """Test inspect prints all attributes when variable_names is None."""
        obj = TriomaClass()
        obj.a = 100
        obj.b = 200

        with patch("sys.stdout", new=StringIO()) as fake_out:
            obj.inspect()
            output = fake_out.getvalue()
            self.assertIn("a: 100", output)
            self.assertIn("b: 200", output)

    def test_inspect_specific_attribute(self):
        """Test inspect with specific variable_names."""
        obj = TriomaClass()
        obj.a = 100
        obj.b = 200

        with patch("sys.stdout", new=StringIO()) as fake_out:
            obj.inspect(variable_names="a")
            output = fake_out.getvalue()
            self.assertIn("a: 100", output)
            self.assertNotIn("b: 200", output)

    def test_inspect_case_insensitive(self):
        """Test inspect matches attribute names case-insensitively."""
        obj = TriomaClass()
        obj.a = 100

        with patch("sys.stdout", new=StringIO()) as fake_out:
            obj.inspect(variable_names="A")
            output = fake_out.getvalue()
            self.assertIn("a: 100", output)

    def test_inspect_nested_trioma_object(self):
        """Test inspect recursively prints nested TriomaClass objects."""
        parent = TriomaClass()
        parent.value = 50
        child = TriomaClass()
        child.nested_attr = 999
        parent.child = child

        with patch("sys.stdout", new=StringIO()) as fake_out:
            parent.inspect()
            output = fake_out.getvalue()
            self.assertIn("value: 50", output)
            self.assertIn("child is a", output)
            self.assertIn("nested_attr: 999", output)


class TestUpdateAttribute(unittest.TestCase):
    """Essential tests for TriomaClass.update_attribute() method."""

    def test_update_direct_attribute(self):
        """Test updating a direct attribute."""
        obj = TriomaClass()
        obj.a = 100
        obj.update_attribute("a", 999)
        self.assertEqual(obj.a, 999)

    def test_update_nested_attribute(self):
        """Test updating attribute in nested TriomaClass object."""
        parent = TriomaClass()
        child = TriomaClass()
        child.a = 100
        parent.child = child

        parent.update_attribute("a", 999)
        self.assertEqual(parent.child.a, 999)

    def test_update_nonexistent_attribute_raises_error(self):
        """Test that updating non-existent attribute raises ValueError."""
        obj = TriomaClass()
        obj.a = 100

        with self.assertRaises(ValueError) as ctx:
            obj.update_attribute("nonexistent", 999)
        self.assertIn("nonexistent", str(ctx.exception))
        self.assertIn("TriomaClass", str(ctx.exception))

    def test_update_n_pipes_propagates_to_nested(self):
        """Test that updating n_pipes propagates to nested TriomaClass objects."""
        parent = TriomaClass()
        parent.n_pipes = 5
        child1 = TriomaClass()
        child1.n_pipes = 5
        child2 = TriomaClass()
        child2.n_pipes = 5
        parent.child1 = child1
        parent.child2 = child2

        parent.update_attribute("n_pipes", 10)

        self.assertEqual(parent.n_pipes, 10)
        self.assertEqual(parent.child1.n_pipes, 10)
        self.assertEqual(parent.child2.n_pipes, 10)

    def test_update_multiple_nesting_levels(self):
        """Test updating attributes through multiple nesting levels."""
        root = TriomaClass()
        middle = TriomaClass()
        leaf = TriomaClass()
        leaf.value = 100
        middle.leaf = leaf
        root.middle = middle

        root.update_attribute("value", 999)
        self.assertEqual(root.middle.leaf.value, 999)


if __name__ == "__main__":
    unittest.main()
