from plotly.subplots import make_subplots
import numpy as np

def plot_subplots_plotly(number_of_rows, number_of_columns, **kwargs):
    fig = make_subplots(rows=number_of_rows, cols=number_of_columns)
    axs = np.array([[(i, j) for j in range(number_of_columns)] for i in range(number_of_rows)])
    return fig, axs