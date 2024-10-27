import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf

class TestToolPrintInColor(unittest.TestCase):
    
    def test_default_color(self):
        """ Test if default color is red """
        try:
            cf.tool_print_in_color("Test string")
        except Exception as e:
            self.fail(f"tool_print_in_color raised an exception with default color: {e}")
    
    def test_valid_color(self):
        """ Test tool_print_in_color with various valid colors """
        valid_colors = ['red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'black']
        for color in valid_colors:
            with self.subTest(color=color):
                try:
                    cf.tool_print_in_color(f"Test string in {color}", color=color)
                except Exception as e:
                    self.fail(f"tool_print_in_color raised an exception with color '{color}': {e}")
    
    def test_invalid_color(self):
        """ Test tool_print_in_color with an invalid color """
        with self.assertRaises(KeyError):
            cf.tool_print_in_color("Test string", color="purple")

    def test_empty_string(self):
        """ Test tool_print_in_color with an empty string """
        try:
            cf.tool_print_in_color("", color="red")
        except Exception as e:
            self.fail(f"tool_print_in_color raised an exception with an empty string: {e}")

if __name__ == '__main__':
    unittest.main()
