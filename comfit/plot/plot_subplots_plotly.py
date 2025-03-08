from typing import Dict, List, Optional, Any, Tuple

import numpy as np
import plotly.graph_objects as go


def plot_subplots_plotly(
    number_of_rows: int,
    number_of_columns: int,
    **kwargs: Any
) -> Tuple[go.Figure, np.ndarray]:
    """Create a plotly figure with subplots.

    Parameters
    ----------
    number_of_rows : int
        Number of rows in the subplot grid.
    number_of_columns : int
        Number of columns in the subplot grid.
    \*\*kwargs : Any
        Additional keyword arguments.

    Returns
    -------
    Tuple[go.Figure, np.ndarray]
        A tuple containing the Plotly figure and axes dictionary.
    """
    fig = go.Figure()
    if number_of_columns == 1:
        axs = np.array([{'row': i+1, 'nrows': number_of_rows, 'col': 1, 'ncols': number_of_columns, 'colorbar': False} for i in range(number_of_rows)])
    elif number_of_rows == 1:
        axs = np.array([{'row': 1, 'nrows': number_of_rows, 'col': i+1, 'ncols': number_of_columns, 'colorbar': False} for i in range(number_of_columns)])
    else:
        axs = np.array([[{'row': i+1, 'nrows': number_of_rows, 'col': j+1, 'ncols': number_of_columns, 'colorbar': False} for j in range(number_of_columns)] for i in range(number_of_rows)])
    return fig, axs