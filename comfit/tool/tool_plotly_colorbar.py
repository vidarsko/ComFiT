from typing import Union, Dict, Any, Literal

# General packages
import numpy as np
import plotly.graph_objects as go

# Local packages
from comfit.tool import tool_generate_numbers_between


superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
                '-': '⁻'}


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

        numbers_between = tool_generate_numbers_between(cmin, cmax)

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
        # ticktext = [tool_format_tick_value(val) for val in tickvals]

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
