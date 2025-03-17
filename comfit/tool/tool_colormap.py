from typing import List, Tuple, Union, Literal

# General packages
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

def tool_colormap(
        colormap_string : str, 
        plot_lib : Literal['plotly', 'matplotlib'] = 'plotly'
        ) -> Union[mcolors.LinearSegmentedColormap, List[Tuple[float, str]]]:
    """Get a colormap for a plot.

    Parameters
    ----------
    colormap_string : str
        The name of the colormap. See options at
        https://comfitlib.com/Plotting/
    plot_lib : str, optional
        The plotting library to use, either 'plotly' or 'matplotlib', by default 'plotly'
    
    Returns
    -------
    Union[mcolors.LinearSegmentedColormap, List[Tuple[float, str]]]
        For plotly: List of tuples with (position, rgb color string)
        For matplotlib: A matplotlib colormap object
    """
    
    
    # Custom colormaps
    sunburst_colors = [  (0.24, 0.15, 0.66),
                    (0.20, 0.48, 0.99),
                    (0.07, 0.74, 0.72),
                    (0.79, 0.75, 0.15),
                    (0.98, 0.98, 0.08)]


    if colormap_string == 'angle':
        # Create hues ranging from cyan to red and back to cyan
        # Note: In the HSV color space, cyan is at 0.5 and red is at 0.
        hues = np.linspace(0.5, 0, 128).tolist() + np.linspace(1, 0.5, 128).tolist()

        # Reverse the order of the hues
        hues = hues[::-1]

        # Convert these hues to RGB colors with full saturation and value
        colors = [mcolors.hsv_to_rgb([h, 1, 1]) for h in hues]

        colormap = mcolors.LinearSegmentedColormap.from_list("custom_hsv", colors)

        if plot_lib == 'plotly':
            # # Convert the colormap to Plotly format
            custom_hsv = []
            for i in range(colormap.N):
                rgb = colormap(i)[:3]
                custom_hsv.append([i / (colormap.N - 1), f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})'])
            return custom_hsv
        else:
            return colormap

    elif colormap_string == 'angle_r':
        # Create hues ranging from cyan to red and back to cyan
        # Note: In the HSV color space, cyan is at 0.5 and red is at 0.
        hues = np.linspace(0.5, 0, 128).tolist() + np.linspace(1, 0.5, 128).tolist()

        # Convert these hues to RGB colors with full saturation and value
        colors = [mcolors.hsv_to_rgb([h, 1, 1]) for h in hues]

        colormap = mcolors.LinearSegmentedColormap.from_list("custom_hsv_r", colors)

        if plot_lib == 'plotly':
            # # Convert the colormap to Plotly format
            custom_hsv = []
            for i in range(colormap.N):
                rgb = colormap(i)[:3]
                custom_hsv.append([i / (colormap.N - 1), f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})'])
            return custom_hsv
        else:
            return colormap

    elif colormap_string == 'bluewhitered':
        # Create hues ranging from cyan to red and back to cyan
        # Define the blue-white-red colormap
        colors = [(0, 'blue'), (0.5, 'white'), (1, 'red')]

        if plot_lib == 'plotly':
            return colors
        else:
            return mcolors.LinearSegmentedColormap.from_list('bluewhitered', colors)
    
    elif colormap_string == 'bluewhitered_r':
        # Create hues ranging from cyan to red and back to cyan
        # Define the blue-white-red colormap
        colors = [(0, 'red'), (0.5, 'white'), (1, 'blue')]

        if plot_lib == 'plotly':
            return colors
        else:
            return mcolors.LinearSegmentedColormap.from_list('bluewhitered_r', colors)

    elif colormap_string == 'sunburst':

        colors = sunburst_colors
        
        if plot_lib == 'plotly':
            cm_data = [(i/4, f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})') 
                      for i, rgb in enumerate(colors)]
            return cm_data

        else:
            positions = [0, 0.25, 0.5, 0.75, 1]
            return mcolors.LinearSegmentedColormap.from_list('sunburst', list(zip(positions, colors)))

    elif colormap_string == 'sunburst_r':

        colors = sunburst_colors[::-1]

        if plot_lib == 'plotly':
            cm_data = [(i/4, f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})') 
                      for i, rgb in enumerate(colors)]
            return cm_data
        elif plot_lib == 'matplotlib':
            positions = [0, 0.25, 0.5, 0.75, 1]
            return mcolors.LinearSegmentedColormap.from_list('sunburst_r', list(zip(positions, colors)))

    elif colormap_string == 'winter':

        if plot_lib == 'matplotlib':
            return plt.get_cmap('winter')
        else:
            # Get the matplotlib colormap
            cmap = plt.get_cmap('winter')
            # Convert to plotly format by sampling N points and converting to RGB strings
            N = 255
            plotly_colors = [(i/(N-1), f'rgb({int(255*r)},{int(255*g)},{int(255*b)})') 
                             for i, (r,g,b) in enumerate([cmap(i)[:3] for i in range(N)])]
            return plotly_colors

    elif colormap_string == 'winter_r':
        
        if plot_lib == 'matplotlib':
            return plt.get_cmap('winter_r')
        else:
            # Get the matplotlib colormap
            cmap = plt.get_cmap('winter_r')
            # Convert to plotly format by sampling N points and converting to RGB strings
            N = 255
            plotly_colors = [(i/(N-1), f'rgb({int(255*r)},{int(255*g)},{int(255*b)})') 
                             for i, (r,g,b) in enumerate([cmap(i)[:3] for i in range(N)])]
            return plotly_colors

    else:
        if plot_lib == 'plotly':
            return colormap_string
        else:
            return plt.get_cmap(colormap_string)