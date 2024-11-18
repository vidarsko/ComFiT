import numpy as np
import matplotlib.pyplot as plt
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib

def plot_field_velocity_and_director_matplotlib(self, field, velocity, director, **kwargs):
    """Plot the fields, velocity, and director field in 2 dimensions

    Args:
        field (ndarray): The field to be plotted.
        velocity (ndarray): The velocity to be plotted.
        director (ndarray): The director to be plotted.
        **kwargs: Keyword arguments for the plot.
            See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.

    Returns:
        A tuple consisting of
            - The figure
            - The axes with the plotted field, velocity, and director. ax (Axes)

    Raises:
        Exception: If the dimension is other than 2.
    """
    if field.dtype == bool:
        field = field.astype(float)

    # Check if the vector field is complex
    if np.iscomplexobj(field):
        print(
                "\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(field)))
    field = np.real(field)

    # Check if an axis object is provided
    fig = kwargs.get('fig', plt.gcf())
    ax = kwargs.get('ax', None)

    # Kewyord arguments
    colorbar = kwargs.get('colorbar', True)

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    if self.dim == 2:

        # Keyword arguments particular to the 2D case
        kwargs['grid'] = kwargs.get('grid', False)

        if ax == None:
            fig.clf()
            ax = plt.gcf().add_subplot(111)

        colormap = kwargs.get('colormap', 'viridis')

        if colormap == 'bluewhitered':
            colormap = tool_colormap_bluewhitered()

        elif colormap == 'sunburst':
            colormap = tool_colormap_sunburst()
        else:
            colormap = plt.get_cmap(colormap)

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')
        pcm = ax.pcolormesh(X / self.a0, Y / self.a0, field, shading='gouraud', cmap=colormap)

        xlim = [self.xmin, self.xmax - self.dx]
        ylim = [self.ymin, self.ymax - self.dy]

        limits_provided = False
        if 'xlim' in kwargs:
            xlim = kwargs['xlim']
            limits_provided = True
        else:
            if 'xmin' in kwargs:
                xlim[0] = kwargs['xmin']
                limits_provided = True

            if 'xmax' in kwargs:
                xlim[1] = kwargs['xmax']
                limits_provided = True

        if 'ylim' in kwargs:
            ylim = kwargs['ylim']
            limits_provided = True
        else:
            if 'ymin' in kwargs:
                ylim[0] = kwargs['ymin']
                limits_provided = True

            if 'ymax' in kwargs:
                ylim[1] = kwargs['ymax']
                limits_provided = True
        if limits_provided:
            region_to_plot = np.zeros(self.dims).astype(bool)
            region_to_plot[(xlim[0] <= X) * (X <= xlim[1]) * (ylim[0] <= Y) * (Y <= ylim[1])] = True
            vlim = [np.min(field[region_to_plot]), np.max(field[region_to_plot])]

        else:
            vlim = [np.min(field), np.max(field)]

        # Set the value limitses
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
        else:
            if 'vmin' in kwargs:
                vlim[0] = kwargs['vmin']
            if 'vmax' in kwargs:
                vlim[1] = kwargs['vmax']

        if vlim[1] - vlim[0] < 1e-10:
            vlim = [vlim[0] - 0.05, vlim[1] + 0.05]

        pcm.set_clim(vmin=vlim[0], vmax=vlim[1])

        if 'vlim_symmetric' in kwargs:
            vlim_symmetric = kwargs['vlim_symmetric']
            if vlim_symmetric:
                cmax = abs(field).max()
                cmin = -cmax
                pcm.set_clim(vmin=cmin, vmax=cmax)

        colorbar = kwargs.get('colorbar', True)

        if colorbar:
            cbar = plt.colorbar(pcm, ax=ax)



        ax.streamplot(X.T, Y.T, (velocity[0]).T, (velocity[1]).T, color='w')
        ax.quiver(X, Y, director[0], director[1], headwidth=0, scale=50)
        ax.quiver(X, Y, -director[0], -director[1], headwidth=0, scale=50)
        ax.set_aspect('equal')

    else:
        raise Exception("This plotting function is currently only implemented in 2D! ")

    kwargs['fig'] = fig
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax