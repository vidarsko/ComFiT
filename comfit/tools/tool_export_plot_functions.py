"""

This module provides tools for exporting plots.
"""


import matplotlib.pyplot as plt
import comfit as cf

from tqdm import tqdm


def tool_export_rotating_plot(ax=None,step=5):


    if ax is None:
        ax = plt.gca()


    for n in tqdm(range(0,720//step),desc="Saving images for animation"):

        ax.view_init(elev=30, azim=-60+n*step)
        cf.tool_save_plot(n)

    cf.tool_make_animation_movie(n,name="rotating_plot.mp4",fps=12)