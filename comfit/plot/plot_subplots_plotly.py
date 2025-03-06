import plotly.graph_objects as go
import numpy as np

def plot_subplots_plotly(number_of_rows, number_of_columns, **kwargs):
    fig = go.Figure()
    if number_of_columns == 1:
        axs = np.array([{'row': i+1, 'nrows': number_of_rows, 'col': 1, 'ncols': number_of_columns, 'colorbar': False} for i in range(number_of_rows)])
    elif number_of_rows == 1:
        axs = np.array([{'row': 1, 'nrows': number_of_rows, 'col': i+1, 'ncols': number_of_columns, 'colorbar': False} for i in range(number_of_columns)])
    else:
        axs = np.array([[{'row': i+1, 'nrows': number_of_rows, 'col': j+1, 'ncols': number_of_columns, 'colorbar': False} for j in range(number_of_columns)] for i in range(number_of_rows)])
    return fig, axs