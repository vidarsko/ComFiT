import matplotlib.pyplot as plt
import imageio
from moviepy.editor import ImageSequenceClip
import os

def tool_save_plot(counter):
    plt.savefig(f'plot_{counter}.png')


def tool_make_animation(counter):
    # List of saved plot filenames
    image_files = [f'plot_{counter}.png' for counter in range(counter+1)]

    # Create the video clip from the image files
    video_clip = ImageSequenceClip(image_files, fps=24)

    # Save the video
    video_clip.write_videofile('output_video.mp4')

    # Delete the png files
    for file in image_files:
        os.remove(file)

