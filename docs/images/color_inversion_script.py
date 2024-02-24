from PIL import Image, ImageOps, ImageSequence
import os

def delete_previous_conversions(folder_path, suffix='-colorinverted'):
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(('.png', '.gif')) and suffix in filename.lower():
            os.remove(os.path.join(folder_path, filename))
            print(f"Deleted {filename}")

def invert_colors(image):
    if image.mode == 'RGBA':
        r,g,b,a = image.split()
        rgb_image = Image.merge('RGB', (r,g,b))
        inverted_image = ImageOps.invert(rgb_image)
        r2,g2,b2 = inverted_image.split()
        return Image.merge('RGBA', (r2,g2,b2,a))
    elif image.mode == 'P':
        image = image.convert('RGB')
        return ImageOps.invert(image)
    else:
        return ImageOps.invert(image)

def invert_image_colors(image_path, output_path):
    with Image.open(image_path) as img:
        if img.format == 'GIF':
            frames = []
            for frame in ImageSequence.Iterator(img):
                inverted_frame = invert_colors(frame)
                frames.append(inverted_frame.convert('P', palette=Image.ADAPTIVE))  # Convert back to P mode with an adaptive palette
            frames[0].save(output_path, save_all=True, append_images=frames[1:], loop=0, duration=img.info.get('duration', 100), optimize=False)
        else:
            inverted_image = invert_colors(img)
            inverted_image.save(output_path)

def process_folder(folder_path):
    delete_previous_conversions(folder_path)
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(('.png', '.gif')):
            base_name, extension = os.path.splitext(filename)
            inverted_filename = f"{base_name}-colorinverted{extension}"
            image_path = os.path.join(folder_path, filename)
            output_path = os.path.join(folder_path, inverted_filename)
            invert_image_colors(image_path, output_path)
            print(f"Processed {filename} -> {inverted_filename}")

# Example usage - replace 'path/to/folder' with the actual folder path
process_folder('./docs/images/')
