import numpy as np
import plotly.graph_objects as go
from comfit.tool import tool_colormap
import math


superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
                '-': '⁻'}

def format_tick_value(val):
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


def generate_numbers_between(cmin, cmax):
    if cmin > cmax:
        raise ValueError("cmin must be less than cmax")

    delta = (cmax - cmin)

    if delta == 0:
        return [0], 0
    
    delta_exp = np.floor(np.log10(max(abs(cmin), abs(cmax))))

    tick_min = np.ceil(cmin/10**delta_exp)*10**delta_exp
    tick_max = np.floor(cmax/10**delta_exp)*10**delta_exp

    number_of_steps = np.floor(delta/10**delta_exp)
    if number_of_steps <= 5:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 10**delta_exp)
    elif number_of_steps <= 10:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 2*10**delta_exp)
    elif number_of_steps <= 15:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 3*10**delta_exp)
    elif number_of_steps <= 20:
        numbers_between = np.arange(tick_min, tick_max+10**delta_exp, 4*10**delta_exp)

    return numbers_between, delta_exp
     
        

def tool_plotly_colorbar(ax, type='normal'):
    """Add a colorbar for the angle field to a plotly figure.

    Args:
        ax: The (row,column) index of the subplot to add the colorbar to.
    
    Returns:
        plotly.graph_objects.Scatter: A scatter plot with no data that can be used as a colorbar.
    """
 
    if type == 'normal':
        

        cmin = ax['vmin']
        cmax = ax['vmax']

        numbers_between = generate_numbers_between(cmin, cmax)

        tickvals, delta_exp = numbers_between
        
        ticktext = [round(tickval/10**delta_exp) for tickval in tickvals]

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