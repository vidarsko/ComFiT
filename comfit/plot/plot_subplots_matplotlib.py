# Standard library imports
from typing import Tuple, Any

# Third-party library imports
import matplotlib.pyplot as plt
import numpy as np


def plot_subplots_matplotlib(
    number_of_rows: int,
    number_of_columns: int,
    **kwargs: Any
    ) -> Tuple[plt.Figure, plt.Axes]:
    """Create a figure with a grid of subplots.

    Parameters
    ----------
    number_of_rows : int
        Number of rows in the subplot grid.
    number_of_columns : int
        Number of columns in the subplot grid.
    \*\*kwargs : Any
        Additional keyword arguments passed to plt.subplots().

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure object and array of Axes objects representing the subplots.
    """
    fig, axs = plt.subplots(number_of_rows, number_of_columns, **kwargs)
    return fig, axs