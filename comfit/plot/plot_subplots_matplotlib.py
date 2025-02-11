import matplotlib.pyplot as plt
import numpy as np

def plot_subplots_matplotlib(number_of_rows, number_of_columns, *args):
    """Makes subplots from a series of matplotlib figures."""
    if len(args) > number_of_rows * number_of_columns:
        raise Exception("Too many subplots for given rows and columns.")

    fig, axes = plt.subplots(number_of_rows, number_of_columns, figsize=(15, 10))

    for i, (old_fig, old_ax) in enumerate(args):
        row = i // number_of_columns
        col = i % number_of_columns

        if isinstance(axes, np.ndarray): # Check if axes is an array (multiple subplots)
            ax = axes[row, col] if number_of_rows > 1 or number_of_columns > 1 else axes[row] # Handle single row/column case
        else:
            ax = axes

        for line in old_ax.lines:
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            ax.plot(xdata, ydata, **line.get_kwargs())  # Recreate the line

        # for patch in old_ax.patches:
        #     ax.add_patch(patch)

        # for image in old_ax.images:
        #     ax.imshow(image.get_data(), extent=image.get_extent(), origin=image.get_origin()) # Example for images

    plt.tight_layout()
    return fig