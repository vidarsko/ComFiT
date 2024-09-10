import matplotlib.pyplot as plt
import numpy as np

import matplotlib 
import mpl_toolkits
import plotly.graph_objects as go

from typing import Union

def tool_set_plot_axis_properties_matplotlib(self, **kwargs) -> Union[matplotlib.axes.Axes, go.Figure]:
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
    
    ####################
    #### MATPLOTLIB ####
    ####################
    if plot_lib == 'matplotlib':
        
        ##### AXES #####
        ax = kwargs.get('ax', plt.gca())

        ##### SIZE #####
        if size is not None:
            print("\033[91mWarning: The size keyword is not valid for matplotlib plots.\033[0m")

        ##### TICKS #####
        if xticks is not None:
            ax.set_xticks(xticks)
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        
        if yticks is not None:
            ax.set_yticks(yticks)
        if yticklabels is not None:
            ax.set_yticklabels(yticklabels)

        if zticks is not None:
            ax.set_zticks(zticks)
        if zticklabels is not None:
            ax.set_zticklabels(zticklabels)

        ##### TITLE #####
        if title is not None:
            ax.set_title(title)

        ##### SUPTITLE #####
        if suptitle is not None:
            ax.get_figure().suptitle(suptitle)

        ##### AXIS LABELS #####
        ax.set_xlabel(xlabel)
        
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        
        if zlabel is not None:
            ax.set_zlabel(zlabel)

        ##### AXIS LIMITS #####
        if isinstance(ax, mpl_toolkits.mplot3d.Axes3D):
            ax.set_xlim3d(xlim[0], xlim[1])

            if ylim is not None:
                ax.set_ylim3d(ylim[0], ylim[1])

            if zlim is not None:
                ax.set_zlim3d(zlim[0], zlim[1])

        else:
            ax.set_xlim(xlim[0], xlim[1])

            if ylim is not None:
                ax.set_ylim(ylim[0], ylim[1])

            if zlim is not None:
                ax.set_zlim(zlim[0], zlim[1])

        ##### GRID #####
        ax.grid(grid)

        ##### AXIS ASPECT RATIO #####
        if axis_equal:
            ax.set_aspect('equal')
    

    ####################
    ###### PLOTLY ######
    ####################
    elif plot_lib == 'plotly':

        ##### FIGURE #####
        fig = kwargs.get('fig', go.Figure())

        if size is None:
            fig.update_layout(width=500, height=500)
        else:
            fig.update_layout(width=size[0], height=size[1])

        ##### TICKS #####
        if xticks is not None:
            fig.update_layout(xaxis=dict(tickvals=xticks))
        if xticklabels is not None:
            fig.update_layout(xaxis=dict(ticktext=xticklabels))

        if yticks is not None:
            fig.update_layout(yaxis=dict(tickvals=yticks))
        if yticklabels is not None:
            fig.update_layout(yaxis=dict(ticktext=yticklabels))

        if zticks is not None:
            fig.update_layout(zaxis=dict(tickvals=zticks))
        if zticklabels is not None:
            fig.update_layout(zaxis=dict(ticktext=zticklabels))

        ##### TITLE #####
        if title is not None:
            fig.update_layout(title_text=suptitle)

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
        fig.update_layout(xaxis_range=xlim)

        if ylim is not None:
            fig.update_layout(yaxis_range=ylim)

        # if zlim is not None:
        #     fig.update_layout(zaxis_range=zlim)

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