
from mpl_toolkits.mplot3d import Axes3D

def tool_matplotlib_define_3D_plot_ax(fig,ax):
    """
    Defines the axes for a 3D plot in matplotlib.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes._subplots.Axes3DSubplot
        Axes object.
    
    Returns
    -------
    ax : matplotlib.axes._subplots.Axes3DSubplot
        Axes object.
    """

    if ax is not None and ax is not isinstance(ax, Axes3D):
        # Get row and column
        subplotspec = ax.get_subplotspec()
        gridspec = subplotspec.get_gridspec()
        row = subplotspec.rowspan.start
        col = subplotspec.colspan.start

        ax.remove()

        fig.add_subplot(gridspec[row, col], projection='3d')
        ax = fig.get_axes()[-1]

    if ax is None:
        ax = fig.add_subplot(111, projection='3d')
    
    return ax