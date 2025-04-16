from typing import List, Optional, Tuple, Union

# Local imports
from .tool_print_in_color import tool_print_in_color

def tool_configure_axis(
        dim: int, 
        name: str, 
        xlim: Optional[List[float]]=None, 
        xmin: Optional[float]=None, 
        xmax: Optional[float]=None, 
        xRes: Optional[int]=None, 
        dx: Optional[float]=None
        ) -> Tuple[List[float], float, float, int, float]:
    """Configure a single axis (x, y, or z) based on the provided parameters.

    When constructing a grid, providing all parameters will over-specify the axis.
    This function takes care of the hierarchy of values and determines the axis limits,
    minimum value, maximum value, resolution, and step size.

    The hierarchy of values is as follows:
    xlim > xmin > xmax > xRes > dx

    Parameters
    ----------
    dim : int
        Dimension of the system (1, 2, or 3)
    name : str
        Name of the axis ('x', 'y', 'z') for error messages
    xlim : list of float, optional
        A 2-element list [min, max] for the axis limits
    xmin : float, optional
        Minimum value of the axis
    xmax : float, optional 
        Maximum value of the axis
    xRes : int, optional
        Resolution (number of points) along the axis
    dx : float, optional
        Step size (spacing between points) along the axis

    Returns
    -------
    tuple
        (xlim, xmin, xmax, xRes, dx) where:
        - xlim is a list of [min, max] values
        - xmin is the minimum value (float)
        - xmax is the maximum value (float) 
        - xRes is the resolution (int)
        - dx is the step size (float)

    Notes
    -----
    Even though the function takes parameters 'x'-specific parameters, 
    it is used for 'y' and 'z' as well.
    """

    # Find the axis limits
    xlim_YUP_provided = xlim is not None
    xlim_NOT_provided = xlim is None
    
    xmin_YUP_provided = xmin is not None
    xmin_NOT_provided = xmin is None

    xmax_YUP_provided = xmax is not None
    xmax_NOT_provided = xmax is None

    xRes_YUP_provided = xRes is not None
    xRes_NOT_provided = xRes is None

    dx_YUP_provided = dx is not None
    dx_NOT_provided = dx is None

    default_xRes = 101
    default_dx = 1.0
    default_xmin = 0
    default_xmax = 101

    if xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. no values provided

        xmin = default_xmin
        xmax = default_xmax 
        xRes = default_xRes
        dx = default_dx

        if (name == 'y' and dim < 2) or (name == 'z' and dim < 3):
            xmax = 1
            xRes = 1
            
        xlim = [xmin, xmax-dx]
        
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. dx = 0.1

        xmin = default_xmin
        xmax = default_xmax
        xRes = round((xmax - xmin) / dx)
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xRes = 56

        xmin = default_xmin
        xmax = xRes
        dx = (xmax - xmin) / xRes
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xRes = 56, dx = 0.1
        
        xmin = default_xmin
        xmax = xmin + xRes * dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xmax = 56
        
        xmin = default_xmin if xmax > 0 else xmax-default_xmax
        xRes = default_xRes
        dx = (xmax - xmin) / xRes
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xmax = 21, dx = 0.1

        xmin = default_xmin if xmax > 0 else xmax-default_xmax
        xRes = round((xmax - xmin) / dx)
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xmax = 21, xRes = 56

        xmin = default_xmin if xmax > 0 else xmax-xRes
        dx = (xmax - xmin) / xRes
        xlim = [xmin, xmax-dx]

    elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xmax = 21, xRes = 56, dx = 0.1

        xmin = xmax - xRes * dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xmin = -5

        xmax = xmin+default_xmax
        xRes = default_xRes
        dx = default_dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xmin = -5, dx = 0.1

        xRes = default_xRes
        xmax = xmin + xRes * dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xmin = -5, xRes = 56

        dx = default_dx
        xmax = xmin + xRes * dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xmin = -5, xRes = 56, dx = 0.1

        xmax = xmin + xRes * dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xmin = -5, xmax = 21

        xRes = xmax - xmin
        dx = default_dx
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xmin = -5, xmax = 21, dx = 0.1

        xRes = round((xmax - xmin) / dx)
        dx = (xmax - xmin) / xRes 
        xlim = [xmin, xmax-dx]

    elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xmin = -5, xmax = 21, xRes = 56

        dx = (xmax - xmin) / xRes
        xlim = [xmin, xmax-dx]
    
    elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xmin = -5, xmax = 21, xRes = 56, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}min, {name}max, and {name}Res trumps providing d{name}.', 'yellow')
        dx = (xmax - xmin) / xRes
        xlim = [xmin, xmax-dx]
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21]

        xRes = default_xRes
        dx = (xlim[1]-xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1]+dx
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], dx = 0.1

        xRes = round((xmax - xmin) / dx)+1
        dx = (xlim[1]-xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1]+dx
    
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xRes = 56
        dx = (xlim[-1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1]+dx
        
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xRes = 56, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim and xRes trumps providing d{name}.', 'yellow')
        
        dx = (xlim[-1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1]+dx
        
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmax = 21

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (default_xRes-1)
        xRes = default_xRes

        xmin = xlim[0]
        xmax = xlim[1]+dx

    elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmax = 21, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

        xRes = round(xlim[1] - xlim[0]) / dx + 1
        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmax = 21, xRes = 56

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmax = 21, xRes = 56, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')
        tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx

    elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmin = -5

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (default_xRes-1)
        xRes = default_xRes

        xmin = xlim[0]
        xmax = xlim[1] + dx

    elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmin = -5, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

        xRes = round(xlim[1] - xlim[0]) / dx + 1
        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx

    elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xRes = 56

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xRes = 56, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')
        tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xmax = 21

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (default_xRes-1)
        xRes = default_xRes

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

        xRes = round(xlim[1] - xlim[0] / dx) + 1
        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, xRes = 56

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
        # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, xRes = 56, dx = 0.1

        tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')
        tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

        dx = (xlim[1] - xlim[0]) / (xRes-1)

        xmin = xlim[0]
        xmax = xlim[1] + dx
    
    return xlim, xmin, xmax, xRes, dx