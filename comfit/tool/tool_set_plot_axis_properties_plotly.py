import numpy as np
import plotly.graph_objects as go

def tool_set_plot_axis_properties_plotly(self, **kwargs):
    """Sets the properties of the axis for a plot.

    Args:
    kwargs: keyword arguments for the axis properties

    Returns:
    The axis object with the properties set.
    """
    # Create dictionaries to store the layout updates
    layout_updates = {}
    
    # For 3D plots, create dictionaries to store the scene updates
    scene_updates = {}
    xaxis_updates = {}
    yaxis_updates = {}
    zaxis_updates = {}

    ##### FIGURE #####
    fig = kwargs.get('fig', go.Figure())

    ##### SIZE #####
    size = kwargs.get('size', None)

    if size is not None:
        layout_updates['width'] = size[0]
        layout_updates['height'] = size[1]

    ##### PLOT NATURE #####
    row = kwargs.get('row', None)
    col = kwargs.get('col', None)

    plot_is_3D = kwargs.get('plot_is_3D', False)

    ##### GRID #####
    grid = kwargs.get('grid', True)
    layout_updates['xaxis_showgrid'] = grid
    layout_updates['yaxis_showgrid'] = grid

    ##### AXIS EQUAL #####
    axis_equal = kwargs.get('axis_equal', False if self.dim==1 else True)

    if axis_equal:
        if plot_is_3D:
            scene_updates['aspectmode'] = 'cube'
        else:
            layout_updates['yaxis'] = dict(scaleanchor='x', scaleratio=1)
    

    ##### TITLE #####
    title = kwargs.get('title', None)

    if title is not None:
        layout_updates['title'] = title

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

    if plot_is_3D:
        xaxis_updates['range'] = xlim
        yaxis_updates['range'] = ylim
        zaxis_updates['range'] = zlim

    else:
        layout_updates['xaxis_range'] = xlim
        layout_updates['yaxis_range'] = ylim

    ##### AXIS LABELS #####
    xlabel = kwargs.get('xlabel', 'x/a₀')
    ylabel = kwargs.get('ylabel', 'y/a₀' if self.dim > 1 else None)
    zlabel = kwargs.get('zlabel', 'z/a₀' if self.dim > 2 else None)

    if plot_is_3D:
        xaxis_updates['title'] = xlabel
        yaxis_updates['title'] = ylabel
        zaxis_updates['title'] = zlabel

    else: #plot is not 3D
        layout_updates['xaxis_title'] = xlabel

        if ylabel is not None:
            layout_updates['yaxis_title'] = ylabel


    ##### TICKS #####
    xticks = kwargs.get('xticks', None)  
    xticklabels = kwargs.get('xticklabels', None)

    yticks = kwargs.get('yticks', None)
    yticklabels = kwargs.get('yticklabels', None)

    zticks = kwargs.get('zticks', None)
    zticklabels = kwargs.get('zticklabels', None)

    
    if xticks is not None:
        xaxis_updates['tickvals'] = xticks
    if xticklabels is not None:
        xaxis_updates['ticktext'] = xticklabels

    if yticks is not None:
        yaxis_updates['tickvals'] = yticks
    if yticklabels is not None:
        yaxis_updates['ticktext'] = yticklabels

    
    if zticks is not None:
        zaxis_updates['tickvals'] = zticks
    if zticklabels is not None:
        zaxis_updates['ticktext'] = zticklabels

    ##### UPDATE SCENE #####
    if plot_is_3D:
        if xaxis_updates:
            scene_updates['xaxis'] = xaxis_updates

        if yaxis_updates:
            scene_updates['yaxis'] = yaxis_updates

        if zaxis_updates:
            scene_updates['zaxis'] = zaxis_updates

    else:
        if xaxis_updates:
            layout_updates['xaxis'] = xaxis_updates

        if yaxis_updates:
            layout_updates['yaxis'] = yaxis_updates


    ##### UPDATE LAYOUT #####
    fig.update_layout(layout_updates)
    if plot_is_3D:
        fig.update_layout(scene = scene_updates)

    
    