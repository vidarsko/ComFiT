import matplotlib.pyplot as plt
import numpy as np

def plot_subplots_matplotlib(number_of_rows, number_of_columns, **kwargs):
    fig, axs = plt.subplots(number_of_rows, number_of_columns)
    return fig, axs