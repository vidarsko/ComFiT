

def plot_save_plotly(counter, fig, **kwargs):
    """
    Saves the current Plotly figure as a PNG file.

    Args:
        counter (int): A unique identifier for the plot image file.
        fig (plotly.graph_objects.Figure): The Plotly figure to save as an image.
        ID (str, optional): A unique identifier for the plot image file. Defaults to None.
    
    Returns: 
        None, saves the plot as a PNG file.
    """

    # Keyword arguments
    ID = kwargs.get('ID', None)

    if ID is None:
        filename = f'plot_{counter}.png'
    else:
        filename = f'plot_{counter}_{ID}.png'

    fig.write_image(filename)