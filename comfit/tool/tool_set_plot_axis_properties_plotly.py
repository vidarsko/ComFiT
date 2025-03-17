from typing import List, Tuple, Union, TYPE_CHECKING, Any
if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

import numpy as np
import plotly.graph_objects as go

def tool_set_plot_axis_properties_plotly(
        self : 'BaseSystem', 
        **kwargs : Any
        ) -> None:
    """Set the properties of the axis for a plot.

    Parameters
    ----------
    kwargs : dict
        Keyword arguments for the axis properties.

    Returns
    -------
    None
        Function updates the layout and scene dictionaries with the axis properties.    
    """
    # Create dictionaries to store the layout updates
    layout_updates = {}
    
    # For 3D plots, create dictionaries to store the scene updates
    scene_updates = {}
    xaxis_updates = {}
    yaxis_updates = {}
    zaxis_updates = {}

    ##### FIGURE #####
    if 'fig' not in kwargs:
        raise ValueError("'fig' parameter is required")
    if 'ax' not in kwargs:
        raise ValueError("'ax' parameter is required")
    
    fig = kwargs['fig']
    ax = kwargs['ax']

    ##### SIZE #####
    size = kwargs.get('size', None)

    if size is not None:
        layout_updates['width'] = size[0]
        layout_updates['height'] = size[1]

    ##### PLOT NATURE #####
    row = ax['row']
    col = ax['col']

    plot_is_3D = kwargs.get('plot_is_3D', False)

    ##### GRID #####
    grid = kwargs.get('grid', True)
    if ax['plot_dimension'] == 2:
        layout_updates[f'{ax["xaxisN"]}_showgrid'] = grid
        layout_updates[f'{ax["yaxisN"]}_showgrid'] = grid

    ##### AXIS EQUAL #####
    axis_equal = kwargs.get('axis_equal', False if self.dim==1 else True)

    if axis_equal:
        if ax['plot_dimension'] == 2:
            layout_updates[ax['yaxisN']] = dict(scaleanchor=ax['xN'], scaleratio=1)
        elif ax['plot_dimension'] == 3:
            scene_updates['aspectmode'] = 'cube'
    

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
    elif 'vlim' in kwargs:
        if self.dim == 1:
            ylim = kwargs['vlim']

    # zlim is specified as a list if dim>2 else as None
    zlim = [self.zmin/self.a0, (self.zmax-self.dz)/self.a0] if self.dim > 2 else None
    if 'zmin' in kwargs:
        zlim[0] = kwargs['zmin'] / self.a0 if self.dim > 2 else kwargs['zmin']
    if 'zmax' in kwargs:
        zlim[1] = kwargs['zmax'] / self.a0 if self.dim > 2 else kwargs['zmax']
    if 'zlim' in kwargs:
        zlim = np.array(kwargs['zlim'])/self.a0 if self.dim > 2 else kwargs['zlim']

    if ax['plot_dimension'] == 2:
        layout_updates[f"{ax['xaxisN']}_range"] = xlim
        layout_updates[f"{ax['xaxisN']}_range"] = ylim
    
    elif ax['plot_dimension'] == 3:
        xaxis_updates['range'] = xlim
        yaxis_updates['range'] = ylim
        zaxis_updates['range'] = zlim
        
    ##### AXIS LABELS #####
    xlabel = kwargs.get('xlabel')
    ylabel = kwargs.get('ylabel', 'y/a₀' if self.dim > 1 else None)
    zlabel = kwargs.get('zlabel', 'z/a₀' if self.dim > 2 else None)

    if ax['plot_dimension'] == 2:
        layout_updates[f"{ax['xaxisN']}_title"] = xlabel
        if ylabel is not None:
            layout_updates[f"{ax['yaxisN']}_title"] = ylabel

    elif ax['plot_dimension'] == 3:
        xaxis_updates['title'] = xlabel
        yaxis_updates['title'] = ylabel
        zaxis_updates['title'] = zlabel
        
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
    if ax['plot_dimension'] == 3:
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
    # if ax is None:

    ##### ADJUST SUBPLOT POSITION #####
    padding = 0.15

    x_domain_start = (ax['col']-1)/ax['ncols']+padding*(ax['col']-1)/ax['ncols']
    x_domain_end = ax['col']/ax['ncols']-padding*(1-ax['col']/ax['ncols'])
    y_domain_start = 1-ax['row']/ax['nrows']+padding*(ax['nrows']-ax['row'])/ax['nrows']
    y_domain_end = 1-(ax['row']-1)/ax['nrows']-padding*(ax['row']-1)/ax['nrows']

    if ax['plot_dimension'] == 2:
        fig.update_layout({ax['xaxisN']: dict( domain=[x_domain_start, 
                                                        x_domain_end],
                                                        anchor=ax['yN'])})

        fig.update_layout({ax['yaxisN']: dict(domain=[y_domain_start, 
                                                      y_domain_end],
                                                            anchor=ax['xN'])})

        fig.update_layout(layout_updates)

    elif ax['plot_dimension'] == 3:
        scene_updates['domain'] = {'x': [x_domain_start, x_domain_end],
                                   'y': [y_domain_start, y_domain_end]}
                                   
        fig.update_layout({ax['sceneN']: scene_updates})

    ##### TITLE #####
    title = kwargs.get('title', None)
    
    if title is not None:
        fig.add_annotation(x=(
                    x_domain_start+0.6*padding*(ax['col']-1)/ax['ncols']
                    +x_domain_end-0.6*padding*(1-ax['col']/ax['ncols']))/2, 
        y= 0.05 + y_domain_end-0.3*padding*(ax['row']-1)/ax['nrows'], 
        xref='paper',
        yref='paper',
        text=title, 
        showarrow=False, 
        align='center')
    

    # else:
    #     dummy_fig = go.Figure()
    #     dummy_fig.update_layout(layout_updates)
    #     if plot_is_3D:
    #         dummy_fig.update_layout(scene=scene_updates)

    #     row = int(ax[0,0])
    #     nrows = int(ax[0,1])
    #     col = int(ax[1,0])
    #     ncols = int(ax[1,1])

    #     if not plot_is_3D:
    #         fig.update_xaxes(dummy_fig.layout.xaxis, row = row, col = col)
    #         fig.update_yaxes(dummy_fig.layout.yaxis, row = row, col = col)
        
    #     else:
    #         fig.update_scenes(dummy_fig.layout.scene, rows = [row], cols = [col])