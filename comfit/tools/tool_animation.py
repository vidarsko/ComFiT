"""

This module provides tools for saving Matplotlib plots as images and creating animations from these images using MoviePy.


Functions:

- tool_save_plot: Saves a Matplotlib plot as an image file.

- tool_make_animation: Creates an animation from a sequence of image files and saves it as a video file.
"""


import matplotlib.pyplot as plt

from moviepy.editor import ImageSequenceClip
import os
from datetime import datetime
import imageio


def tool_save_plot(counter):
    """

    Saves the current Matplotlib plot as a PNG file.


    Input:

    - counter (int): A unique identifier for the plot image file.
    """

    plt.savefig(f'plot_{counter}.png')


def tool_make_animation_movie(counter, name=datetime.now().strftime("%y%m%d_%H%M") + ' - output_video.mp4', fps=24):
    """

    Creates an animation from a series of plot images and saves it as an MP4 video file.


    Input:

    - counter (int): The number of plot images to include in the animation.

    - name (str, optional): The filename for the output video. Defaults to today's date followed by ' - output_video.mp4'.

    - fps (int, optional): The frames per second for the video. Defaults to 24.
    """


    # List of saved plot filenames

    image_files = [f'plot_{counter}.png' for counter in range(counter+1)]


    # Create the video clip from the image files

    video_clip = ImageSequenceClip(image_files, fps=fps)


    # Save the video

    video_clip.write_videofile(name)


    # Delete the png files

    for file in image_files:

        os.remove(file)


def tool_make_animation_gif(counter, name=datetime.now().strftime("%y%m%d_%H%M") + ' - output_animation.gif', fps=24):
    """
    Creates an animation from a series of plot images and saves it as a GIF file.

    Input:
        - counter (int): The number of plot images to include in the animation.
        - name (str, optional): The filename for the output video. Defaults to today's date followed by ' - output_video.mp4'.
        - fps (int, optional): The frames per second for the video. Defaults to 24.
    """

    images = [f'plot_{counter}.png' for counter in range(counter+1)]
    imageio.mimsave(name, images, fps=fps) 
