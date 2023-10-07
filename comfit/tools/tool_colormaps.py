import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors


def tool_colormap_angle():
    # Create hues ranging from cyan to red and back to cyan
    # Note: In the HSV color space, cyan is at 0.5 and red is at 0.
    hues = np.linspace(0.5, 0, 128).tolist() + np.linspace(1, 0.5, 128).tolist()

    # Convert these hues to RGB colors with full saturation and value
    colors = [mcolors.hsv_to_rgb([h, 1, 1]) for h in hues]

    return LinearSegmentedColormap.from_list("custom_hsv", colors)


def tool_colormap_bluewhitered():
    # Create hues ranging from cyan to red and back to cyan
    # Define the blue-white-red colormap
    colors = [(0, 'blue'), (0.5, 'white'), (1, 'red')]
    return mcolors.LinearSegmentedColormap.from_list('blue_white_red', colors)
