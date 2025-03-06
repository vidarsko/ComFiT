
import matplotlib.pyplot as plt

def tool_matplotlib_define_2D_plot_ax(fig, ax):
    """
    Defines the axes for a 2D plot in matplotlib.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes object.
    
    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes object.
    """

    if ax is None:
        ax = fig.add_subplot(111)
    
    return ax