"""
This module provides tools for saving Matplotlib plots as img and creating animations from these img using MoviePy.


Functions:

- tool_save_plot: Saves a Matplotlib plot as an image file.
- tool_make_animation: Creates an animation from a sequence of image files and saves it as a video file.
"""


import matplotlib.pyplot as plt

from moviepy.editor import ImageSequenceClip
import os
from datetime import datetime
from PIL import Image


def tool_save_plot(counter, image_size_inches=(6,5), dpi=100):
    """

    Saves the current Matplotlib plot as a PNG file.

    Args:
        counter (int): A unique identifier for the plot image file.
        image_size (tuple, optional): The size of the image in inches. Defaults to (5,5).
        dpi (int, optional): The resolution of the image in dots per inch. Defaults to 100.
    
    Returns: 
        None, saves the plot as a PNG file.
    """
    plt.gcf().set_size_inches(image_size_inches)
    plt.savefig(f'plot_{counter}.png',dpi=dpi)


def tool_make_animation_movie(counter, name=None, fps=24):
    """

    Creates an animation from a series of plot img and saves it as an MP4 video file.


    Args:
        counter (int): The number of plot img to include in the animation.
        name (str, optional): The filename for the output video. Defaults to today's date followed by ' - output_video.mp4'.
        fps (int, optional): The frames per second for the video. Defaults to 24.
    
    Returns:
        None, saves the video file.
    """
    
    if name is None:
        name = datetime.now().strftime("%y%m%d_%H%M") + ' - output_video.mp4'
    else:
        name = datetime.now().strftime("%y%m%d_%H%M") + ' - ' + name + '.mp4'

    # List of saved plot filenames
    image_files = [f'plot_{counter}.png' for counter in range(counter+1)]

    # Create the video clip from the image files
    video_clip = ImageSequenceClip(image_files, fps=fps)

    # Save the video
    video_clip.write_videofile(name)

    # Delete the png files
    for file in image_files:
        os.remove(file)


def tool_make_animation_gif(counter, name=None, fps=24):
    """
    Creates an animation from a series of plot img and saves it as a GIF file.

    Args:
        counter (int): The number of plot img to include in the animation.
        name (str, optional): The filename for the output video. Defaults to today's date followed by ' - output_video.mp4'.
        fps (int, optional): The frames per second for the video. Defaults to 24.

    Returns:
        None, saves the GIF file.
    """

    if name is None:
        name = datetime.now().strftime("%y%m%d_%H%M") + ' - output_animation.gif'
    else:
        name = datetime.now().strftime("%y%m%d_%H%M") + ' - ' + name + '.gif'

    image_files = [f'plot_{counter}.png' for counter in range(counter+1)]
    images = [Image.open(image_file) for image_file in image_files]

    # Find the smallest image dimensions
    min_width = min(image.width for image in images)
    min_height = min(image.height for image in images)

    # Crop images if necessary and issue a warning
    cropped_images = []
    for image in images:
        if image.width != min_width or image.height != min_height:
            print(f"Warning: Image sizes inconsistent. Cropping image {image.filename} to {min_width}x{min_height}. Consider reviewing how plots are saved mid loop. ")
            image = image.crop((0, 0, min_width, min_height))
        cropped_images.append(image)

    # Save the images as a GIF
    cropped_images[0].save(name, save_all=True, append_images=cropped_images[1:], duration=1000/fps, loop=0)

    # Delete the png files
    for file in image_files:
        os.remove(file)

