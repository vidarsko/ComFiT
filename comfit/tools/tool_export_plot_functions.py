"""
This module provides tools for exporting plots.
"""
import matplotlib.pyplot as plt
import comfit as cf

from tqdm import tqdm

def tool_export_rotating_plot(ax=None,step=5):
    if ax is None:
        ax = plt.gca()

    for n in tqdm(range(0,360//step),desc="Saving img for animation"):

        ax.view_init(elev=30, azim=-60+n*step)
        cf.tool_save_plot(n)

    cf.tool_make_animation_gif(n,name="rotating_plot",fps=12)