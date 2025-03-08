from typing import Optional

def tool_print_in_color(
        string : str, 
        color : Optional[str] = 'red'):
    """Print a string in color.

    Parameters
    ----------
    string : str
        The string to be printed.
    color : str
        The color of the string. Options are 'red' (default), 'green', 
        'yellow', 'blue', 'magenta', 'cyan', 'white', 'black'.

    Returns
    -------
    None
        Prints the string
    """
    colors = {'red': 31, 'green': 32, 'yellow': 33, 'blue': 34, 'magenta': 35, 'cyan': 36, 'white': 37, 'black': 30}
    print('\033[1;{}m{}\033[0m'.format(colors[color], string))
