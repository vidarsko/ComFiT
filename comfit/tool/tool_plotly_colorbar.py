from typing import Tuple, Union, Dict, Any, Literal

# General packages
import numpy as np
import plotly.graph_objects as go
import math

# Local packages
from comfit.tool import tool_colormap


superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
                '-': '⁻'}

def format_tick_value(
    val: float
    ) -> str:
    """Format a numeric value with appropriate SI prefix and significant figures.
    This function formats numeric values using SI prefixes and adjusts decimal places
    based on the magnitude of the number. For values that don't align with standard
    SI prefixes, scientific notation is used as a fallback.

    Parameters
    ----------
    val : float
        The numeric value to format.

    Returns
    -------
    str
        A string representation of the value with appropriate SI prefix and formatting.
    """
    if val == 0:
        return '0'

    # prefixes = {exp: f'×10{"".join(superscripts[digit] for digit in str(exp))}' for exp in range(-18, 19, 3)}
    # prefix_letters = {
    #      -15: 'f', -12: 'p', -9: 'n', -6: 'µ', -3: 'm',
    #     0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T'
    # }
    # for key in prefix_letters.keys():
    #     prefixes[key] = prefix_letters[key]


    abs_val = abs(val)
    exp = int(np.floor(np.log10(abs_val) / 3) * 3)  # Round to nearest power of 1000
    
    prefixes = {}

    if exp in prefixes:
        scaled_val = val / 10**exp
        # Format based on the magnitude of scaled_val
        if abs(scaled_val) < 10:
            return f'{scaled_val:.2f}{prefixes[exp]}'  # 2.00m
        elif abs(scaled_val) < 100:
            return f'{scaled_val:.1f}{prefixes[exp]}'  # 20.0m
        else:
            return f'{scaled_val:.0f}{prefixes[exp]}'  # 200m
    else:
        return f'{val:.0e}'  # Fallback to scientific notation



def generate_numbers_between(
        cmin: float, 
        cmax: float
        ) -> Tuple[np.ndarray, float]:
    """Generate evenly spaced numbers between given minimum and maximum values.

    This function generates a sequence of numbers between cmin and cmax with
    appropriate step sizes based on the range magnitude. The step size is 
    adjusted to maintain a readable number of ticks.

    Parameters
    ----------
    cmin : float
        Minimum value of the range.
    cmax : float
        Maximum value of the range.

    Returns
    -------
    tuple[np.ndarray, float]
        A tuple containing:
        - np.ndarray: Array of evenly spaced numbers between cmin and cmax
        - float: The exponent of the range magnitude (delta_exp)

    Raises
    ------
    ValueError
        If cmin is greater than cmax.

    Notes
    -----
    The step size is determined based on the range magnitude:
    - If range <= 5 units: step = 10^delta_exp
    - If range <= 10 units: step = 2×10^delta_exp
    - If range <= 15 units: step = 3×10^delta_exp
    - If range <= 20 units: step = 4×10^delta_exp
    """
    if cmin > cmax:
        raise ValueError("cmin must be less than cmax")

    delta = (cmax - cmin)

    if delta == 0:
        return [0], 0, False
    
    delta_exp = np.floor(np.log10(max(abs(cmin), abs(cmax))))

    tick_min = np.ceil(cmin/10**delta_exp)*10**delta_exp
    tick_max = np.floor(cmax/10**delta_exp)*10**delta_exp
    
    number_of_steps = np.floor(delta/10**delta_exp)
    
    reiterated = False
    if number_of_steps < 1:
        cmid = (cmin + cmax) / 2
        numbers_between, delta_exp, reiterated = generate_numbers_between(cmin-cmid, cmax-cmid)
        reiterated = True

    elif number_of_steps <= 5:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 10**delta_exp)
    elif number_of_steps <= 10:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 2*10**delta_exp)
    elif number_of_steps <= 15:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 3*10**delta_exp)
    elif number_of_steps <= 20:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 4*10**delta_exp)

    return numbers_between, delta_exp, reiterated
     
        

def tool_plotly_colorbar(
        ax : Dict[str, Any], 
        type : Literal['normal', 'angle'] = 'normal'
        ) -> Union[go.Scatter, go.Scatter3d]:
    """Add a colorbar to a plotly figure subplot with customizable type and placement.

    Parameters
    ----------
    ax : dict
        Dictionary containing subplot information with keys:
        - 'vmin': Minimum value for normal colorbar
        - 'vmax': Maximum value for normal colorbar
        - 'plot_dimension': Integer (2 or 3) indicating plot dimensions
        - 'colormap_object': Plotly colormap object
        - 'row': Row index of subplot
        - 'col': Column index of subplot
        - 'nrows': Total number of rows in figure
        - 'ncols': Total number of columns in figure

    type : str, optional
        Type of colorbar to create, by default 'normal'
        - 'normal': Linear scale with auto-formatted tick values
        - 'angle': Angular scale from -π to π with radian labels

    Returns
    -------
    Union[plotly.graph_objects.Scatter, plotly.graph_objects.Scatter3d]
        A scatter trace containing the colorbar configuration. The trace has no data points
        and is used solely for displaying the colorbar.

    Notes
    -----
    The colorbar placement is automatically calculated based on the subplot position.
    For normal type, tick values are auto-formatted with scientific notation.
    For angle type, tick values are displayed in π radians.
    """
 
    if type == 'normal':
        

        cmin = ax['vmin']
        cmax = ax['vmax']

        numbers_between = generate_numbers_between(cmin, cmax)

        tickvals, delta_exp, reiterated = numbers_between
        ticktext = [round(tickval/10**delta_exp) for tickval in tickvals]

        if reiterated:
            cmid = (cmin + cmax) / 2
            tickvals = [val + cmid for val in tickvals]

        
        if reiterated:
            title = f'{cmid:.1e}' + '+10' + ''.join([superscripts[digit] for digit in str(int(delta_exp))])  
        else:
            title='×10' + ''.join([superscripts[digit] for digit in str(int(delta_exp))])
        
        # tickvals = np.linspace(cmin, cmax, 7)
        # ticktext = [format_tick_value(val) for val in tickvals]

    elif type == 'angle':
        title=None
        cmin = -np.pi
        cmax = np.pi
        tickvals = [-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi]
        ticktext = ['-π', '-2π/3', '-π/3', '0', 'π/3', '2π/3', 'π']
    
    if ax['plot_dimension'] == 2:
        trace = go.Scatter(
                x=[None], y=[None], mode='markers',
                showlegend=False,
                marker=dict(
                    colorscale=ax['colormap_object'],
                    cmin=cmin, 
                    cmax=cmax,
                    colorbar=dict(
                        tickvals=tickvals,
                        ticktext=ticktext
                    )
                ),
                hoverinfo='none'
                )
                
    elif ax['plot_dimension'] == 3:
        trace = go.Scatter3d(
                x=[None], y=[None], z=[None], mode='markers',
                showlegend=False,
                marker=dict(
                    colorscale=ax['colormap_object'],
                    cmin=cmin, 
                    cmax=cmax,
                    colorbar=dict(
                    tickvals=tickvals,
                    ticktext=ticktext
                    )
                ),
                hoverinfo='none'
                )
    
    # The idea here is to place the colorbar in the middle of the subplot
    # It is not perfect yet (Vidar, 13.02.25)
    # padding=0.0

    y_axis_placement = 1-1/ax['nrows']/2 -(ax['row']-1)/ax['nrows'] 
    #0.5 if nrows = 1, 0.33, 0.66 if nrows = 2, 0.25, 0.5, 0.75 if nrows = 3 etc. 
    y_axis_correction = 0 if ax['nrows'] == 1 else -(ax['row']-(1+ax['nrows'])/2)/(ax['nrows']-(1+ax['nrows'])/2) 
    #Goes from 1 (row=1) to -1 (row=nrows) (zero if nrows = 1)
    y_axis_constant = 0
                    
    trace.marker.colorbar.update(
        title=title,
            len=1/ax['nrows'],
            thickness=10,
            xanchor='left',
            x=ax['col']/ax['ncols']-0.1*(1-ax['col']/ax['ncols']),
            yanchor='middle',
            y= y_axis_placement + 0.03*y_axis_correction + y_axis_constant,
        )

    return trace