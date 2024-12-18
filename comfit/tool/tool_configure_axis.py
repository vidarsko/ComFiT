from comfit.tool.tool_print_in_color import tool_print_in_color

def tool_configure_axis(dim, name, xlim=None, xmin=None, xmax=None, xRes=None, dx=None):
        """
        Configure a single axis (x, y, or z) based on the provided parameters.
        Parameters:
            dim (int): Dimension of the system (1, 2, or 3).
            name (str): Name of the axis (e.g., 'x', 'y', 'z') for error messages.
            xlim (list): A 2-element list [min, max] for the axis limits.
            xmin (float): Minimum value of the axis.
            xmax (float): Maximum value of the axis.
            xRes (int): Resolution (number of points) along the axis.
            dx (float): Step size (spacing between points) along the axis.
        Returns:
            xlim, xmin, xmax, xRes, dx
        """

        #Hierarcy of values: lim > min_val > max_val > res > delta

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

        if xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. no values provided

            xmin = 0
            xmax = 101 
            xRes = 101
            dx = 1.0

            if (name == 'y' and dim < 2) or (name == 'z' and dim < 3):
                xmax = 1
                xRes = 1
                
            xlim = [xmin, xmax]
            
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. dx = 0.1

            xmin = 0
            xmax = 101
            xRes = round((xmax - xmin) / dx)
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xRes = 56

            xmin = 0
            xmax = xRes
            dx = (xmax - xmin) / xRes
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xRes = 56, dx = 0.1
            
            xmin = 0
            xmax = xmin + xRes * dx
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xmax = 56
            
            xmin = 0 if xmax > 0 else xmax-101
            xRes = 101
            dx = (xmax - xmin) / xRes
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xmax = 21, dx = 0.1

            xmin = 0 if xmax > 0 else xmax-101
            xRes = round((xmax - xmin) / dx)
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xmax = 21, xRes = 56

            xmin = 0 if xmax > 0 else xmax-xRes
            dx = (xmax - xmin) / xRes
            xlim = [xmin, xmax]

        elif xlim_NOT_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xmax = 21, xRes = 56, dx = 0.1

            xmin = xmax - xRes * dx
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xmin = -5

            xmax = xmin+101
            xRes = 101
            dx = 1.0
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xmin = -5, dx = 0.1

            xRes = 101
            xmax = xmin + xRes * dx
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xmin = -5, xRes = 56

            dx = 1.0
            xmax = xmin + xRes * dx
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xmin = -5, xRes = 56, dx = 0.1

            xmax = xmin + xRes * dx
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xmin = -5, xmax = 21

            xRes = xmax - xmin
            dx = 1.0
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xmin = -5, xmax = 21, dx = 0.1

            xRes = round((xmax - xmin) / dx)
            dx = (xmax - xmin) / xRes 
            xlim = [xmin, xmax]

        elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xmin = -5, xmax = 21, xRes = 56

            dx = (xmax - xmin) / xRes
            xlim = [xmin, xmax]
        
        elif xlim_NOT_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xmin = -5, xmax = 21, xRes = 56, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}min, {name}max, and {name}Res trumps providing d{name}.', 'yellow')
            dx = (xmax - xmin) / xRes
            xlim = [xmin, xmax]
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21]

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = 101
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], dx = 0.1

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = round((xmax - xmin) / dx)
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xRes = 56

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xRes = 56, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim and xRes trumps providing d{name}.', 'yellow')
            
            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmax = 21

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = 101
            dx = (xmax - xmin) / xRes

        elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmax = 21, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = round((xmax - xmin) / dx)
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmax = 21, xRes = 56

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_NOT_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmax = 21, xRes = 56, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}max.', 'yellow')
            tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes

        elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmin = -5

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = 101
            dx = (xmax - xmin) / xRes

        elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmin = -5, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = round((xmax - xmin) / dx)
            dx = (xmax - xmin) / xRes

        elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xRes = 56

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_YUP_provided and xmax_NOT_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xRes = 56, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min.', 'yellow')
            tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xmax = 21

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = 101
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_NOT_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            xRes = round((xmax - xmin) / dx)
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_NOT_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, xRes = 56

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        elif xlim_YUP_provided and xmin_YUP_provided and xmax_YUP_provided and xRes_YUP_provided and dx_YUP_provided:
            # E.g. xlim = [-5, 21], xmin = -5, xmax = 21, xRes = 56, dx = 0.1

            tool_print_in_color(f'Warning: Providing {name}lim trumps providing {name}min and {name}max.', 'yellow')
            tool_print_in_color(f'Warning: Providing {name}Res and {name}lim trumps providing d{name}.', 'yellow')

            xmin = xlim[0]
            xmax = xlim[1]
            dx = (xmax - xmin) / xRes
        
        return xlim, xmin, xmax, xRes, dx