from typing import Tuple, Optional

#General packages
import matplotlib.pyplot as plt
try:
    from moviepy import ImageSequenceClip
except ImportError:
    from moviepy.editor import ImageSequenceClip
import os
from datetime import datetime
from PIL import Image


def tool_make_animation_movie(
        counter : int, 
        name : Optional[str] = None, 
        fps : Optional[int] = 24,
        ID : Optional[str] = None,
        delete_png: Optional[bool] = True
        ) -> None:
    """Make an animation from a series of plot img and save it as an MP4 video file.

    Parameters
    ----------
    counter : int
        The number of plot img to include in the animation.
    name : Optional[str], optional
        The filename for the output video, by default None
    fps : Optional[int], optional
        The frames per second for the video, by default 24
    ID : Optional[str], optional
        A unique identifier for the plot image files, by default None
    delete_png : Optional[bool], optional
        Whether to delete the png files after creating the video, by default True

    Returns
    -------
    None
        Saves the video file.
    """

    if name is None:
        name = datetime.now().strftime("%y%m%d_%H%M") + (' - ' + ID if ID else '') + ' - output_video.mp4'
    else:
        name = datetime.now().strftime("%y%m%d_%H%M") + (' - ' + ID if ID else '') + ' - ' + name + '.mp4'

    # List of saved plot filenames
    if ID is None:
        image_files = [f'plot_{counter}.png' for counter in range(counter+1)]
    else:
        image_files = [f'plot_{counter}_{ID}.png' for counter in range(counter+1)]

    # Create the video clip from the image files
    video_clip = ImageSequenceClip(image_files, fps=fps)

    # Save the video
    video_clip.write_videofile(name)

    # Delete the png files
    if delete_png:
        for file in image_files:
            os.remove(file)


def tool_make_animation_gif(
        counter: int,
        name: Optional[str] = None,
        fps: Optional[int] = 24,
        ID: Optional[str] = None,
        delete_png: Optional[bool] = True
        ) -> None:
    """Make an animation from a series of plot img and save it as a GIF file.

    Parameters
    ----------
    counter : int
        The number of plot img to include in the animation.
    name : Optional[str], optional
        The filename for the output GIF, by default None
    fps : Optional[int], optional
        The frames per second for the GIF, by default 24
    ID : Optional[str], optional
        A unique identifier for the plot image files, by default None
    delete_png : Optional[bool], optional
        Whether to delete the png files after creating the GIF, by default True

    Returns
    -------
    None
        Saves the GIF file.
    """

    if name is None:
        name = datetime.now().strftime("%y%m%d_%H%M%S") + (' - ' + ID if ID else '') + ' - output_animation.gif'
    else:
        name = datetime.now().strftime("%y%m%d_%H%M%S") + (' - ' + ID if ID else '') + ' - ' + name + '.gif'

    if ID is None:
        image_files = [f'plot_{counter}.png' for counter in range(counter+1)]
    else:
        image_files = [f'plot_{counter}_{ID}.png' for counter in range(counter+1)]
        
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
    if delete_png:
        for file in image_files:
            os.remove(file)

