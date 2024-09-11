import matplotlib.pyplot as plt
import numpy as np

import matplotlib 
import mpl_toolkits
import plotly.graph_objects as go

from typing import Union

def tool_set_plot_axis_properties_plotly(self, **kwargs):
    """Sets the properties of the axis for a plot.

    Args:
    kwargs: keyword arguments for the axis properties

    Returns:
    The axis object with the properties set.
    """

    ##### AXIS LIMITS #####
    # xlim is specified as a list
    xlim = [self.xmin/self.a0, (self.xmax-self.dx)/self.a0]
    if 'xmin' in kwargs:
        xlim[0] = kwargs['xmin'] / self.a0
    if 'xmax' in kwargs:
        xlim[1] = kwargs['xmax'] / self.a0
    if 'xlim' in kwargs:
        xlim = np.array(kwargs['xlim']) / self.a0

    # ylim is specified as a list if dim>1 else as None
    ylim = [self.ymin/self.a0, (self.ymax-self.dy)/self.a0] if self.dim > 1 else None
    if 'ymin' in kwargs:
        ylim[0] = kwargs['ymin'] / self.a0 if self.dim > 1 else kwargs['ymin']
    if 'ymax' in kwargs:
        ylim[1] = kwargs['ymax'] / self.a0 if self.dim > 1 else kwargs['ymax']
    if 'ylim' in kwargs:
        ylim = np.array(kwargs['ylim'])/self.a0 if self.dim > 1 else kwargs['ylim']

    # zlim is specified as a list if dim>2 else as None
    zlim = [self.zmin/self.a0, (self.zmax-self.dz)/self.a0] if self.dim > 2 else None
    if 'zmin' in kwargs:
        zlim[0] = kwargs['zmin'] / self.a0 if self.dim > 2 else kwargs['zmin']
    if 'zmax' in kwargs:
        zlim[1] = kwargs['zmax'] / self.a0 if self.dim > 2 else kwargs['zmax']
    if 'zlim' in kwargs:
        zlim = np.array(kwargs['zlim'])/self.a0 if self.dim > 2 else kwargs['zlim']

    ##### GRID AND TITLE #####
    grid = kwargs.get('grid', True)
    axis_equal = kwargs.get('axis_equal', True)
    title = kwargs.get('title', None)
    suptitle = kwargs.get('suptitle', None)

    ##### SIZE #####
    size = kwargs.get('size', None)

    ##### TICKS #####
    xticks = kwargs.get('xticks', None)  
    xticklabels = kwargs.get('xticklabels', None)

    yticks = kwargs.get('yticks', None)
    yticklabels = kwargs.get('yticklabels', None)

    zticks = kwargs.get('zticks', None)
    zticklabels = kwargs.get('zticklabels', None)

    ##### LABELS #####
    xlabel = kwargs.get('xlabel', 'x/a₀')
    ylabel = kwargs.get('ylabel', 'y/a₀' if self.dim > 1 else None)
    zlabel = kwargs.get('zlabel', 'z/a₀' if self.dim > 2 else None)

    ##### PLOT NATURE #####
    plot_is_3D = kwargs.get('plot_is_3D', False)

    ##### PLOT LIBRARY #####
    plot_lib = kwargs.get('plot_lib', self.plot_lib)


    ##### FIGURE #####
    fig = kwargs.get('fig', go.Figure())

    if size is None:
        fig.update_layout(width=500, height=500)
    else:
        fig.update_layout(width=size[0], height=size[1])

    ##### TICKS #####
    if plot_is_3D:
        if xticks is not None:
            fig.update_layout(scene=dict(xaxis=dict(tickvals=xticks)))
        if xticklabels is not None:
            fig.update_layout(scene=dict(xaxis=dict(ticktext=xticklabels)))

        if yticks is not None:
            fig.update_layout(scene=dict(yaxis=dict(tickvals=yticks)))
        if yticklabels is not None:
            fig.update_layout(scene=dict(yaxis=dict(ticktext=yticklabels)))

        if zticks is not None:
            fig.update_layout(scene=dict(zaxis=dict(tickvals=zticks)))
        if zticklabels is not None:
            fig.update_layout(scene=dict(zaxis=dict(ticktext=zticklabels)))

    else:
        if xticks is not None:
            fig.update_layout(xaxis=dict(tickvals=xticks))
        if xticklabels is not None:
            fig.update_layout(xaxis=dict(ticktext=xticklabels))

        if yticks is not None:
            fig.update_layout(yaxis=dict(tickvals=yticks))
        if yticklabels is not None:
            fig.update_layout(yaxis=dict(ticktext=yticklabels))

    ##### TITLE #####
    if title is not None:
        fig.update_layout(title_text=title)

    ##### SUPTITLE #####
    if suptitle is not None:
        print("\033[91mWarning: The suptitle keyword is not valid for plotly plots.\033[0m")

    ##### AXIS LABELS #####
    # Figure is a 3D plot
    if plot_is_3D:
        fig.update_layout(
        scene=dict(
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            zaxis_title=zlabel
        )
    )

    # Figure is not 3D plot
    else:
        fig.update_layout(xaxis_title=xlabel)

    if ylabel is not None:
        fig.update_layout(yaxis_title=ylabel)


    ##### AXIS LIMITS #####
    if plot_is_3D:
        fig.update_layout(scene=dict(
            xaxis=dict(range=xlim),  # Set x-axis range
            yaxis=dict(range=ylim),  # Set y-axis range
            zaxis=dict(range=zlim)   # Set z-axis range
        ))
    else:
        fig.update_layout(xaxis=dict(range=xlim), 
                          yaxis=dict(range=ylim))


    ##### GRID #####

    fig.update_layout(
    xaxis=dict(showgrid=grid),  # Show grid on x-axis
    yaxis=dict(showgrid=grid)   # Show grid on y-axis
    )

    ##### AXIS ASPECT RATIO #####
    if axis_equal:
        fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1)