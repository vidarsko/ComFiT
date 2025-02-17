import numpy as np
import plotly.graph_objects as go
from comfit.tool import tool_colormap_angle

def format_tick_value(val):
    if val == 0:
        return '0'
        
    

    superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                   '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
                   '-': '⁻'}

    prefixes = {exp: f'×10{"".join(superscripts[digit] for digit in str(exp))}' for exp in range(-18, 19, 3)}
    prefix_letters = {
         -15: 'f', -12: 'p', -9: 'n', -6: 'µ', -3: 'm',
        0: '', 3: 'k', 6: 'M', 9: 'G', 12: 'T'
    }
    for key in prefix_letters.keys():
        prefixes[key] = prefix_letters[key]


    
    abs_val = abs(val)
    exp = int(np.floor(np.log10(abs_val) / 3) * 3)  # Round to nearest power of 1000
    
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
        return f'{val:.2e}'  # Fallback to scientific notation


def tool_plotly_colorbar(ax, type='normal'):
    """Add a colorbar for the angle field to a plotly figure.

    Args:
        ax: The (row,column) index of the subplot to add the colorbar to.
    
    Returns:
        plotly.graph_objects.Scatter: A scatter plot with no data that can be used as a colorbar.
    """

    if type == 'normal':
        colormap = 'Viridis'
        cmin = ax['vmin']
        cmax = ax['vmax']
        tickvals = np.linspace(cmin, cmax, 7)
        ticktext = [format_tick_value(val) for val in tickvals]

    elif type == 'angle':
        colormap = tool_colormap_angle(pyplot=True)
        cmin = -np.pi
        cmax = np.pi
        tickvals = [-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi]
        ticktext = ['-π', '-2π/3', '-π/3', '0', 'π/3', '2π/3', 'π']
    

    if ax['plot_dimension'] == 2:
        trace = go.Scatter(
                x=[None], y=[None], mode='markers',
                showlegend=False,
                marker=dict(
                    colorscale=colormap,
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
                    colorscale=colormap,
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
    trace.marker.colorbar.update(
            len=1/ax['nrows'],
            xanchor='right',
            x=ax['col']/ax['ncols'],
            yanchor='middle',
            y=1-1/ax['nrows']/2-(ax['row']-1)/ax['nrows']
        )

    return trace