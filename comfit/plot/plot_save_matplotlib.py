

def plot_save_matplotlib(counter, fig, **kwargs):
    """

    Saves the current Matplotlib plot as a PNG file.

    Args:
        fig (matplotlib.figure.Figure, ax): The Matplotlib figure and axis to save as an image.
        counter (int): A unique identifier for the plot image file.

    Keyword Args:
        ID (str, optional): A unique identifier for the plot image file. Defaults to None.
        image_size (tuple, optional): The size of the image in inches. Defaults to (5,5).
        dpi (int, optional): The resolution of the image in dots per inch. Defaults to 100.
        
    
    Returns: 
        None, saves the plot as a PNG file.
    """

    # Keyword arguments
    ID=kwargs.get('ID', None)
    image_size_inches=kwargs.get('image_size_inches', (6,5))
    dpi=kwargs.get('dpi', 100)

    if ID is None:
        filename = f'plot_{counter}.png'
    else:
        filename = f'plot_{counter}_{ID}.png'

    fig[0].set_size_inches(image_size_inches)
    fig[0].savefig(filename,dpi=dpi)