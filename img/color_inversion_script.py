from PIL import Image, ImageOps, ImageSequence, ImageChops
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

def change_black_to_color(image, new_color):
    data = image.getdata()
    
    newData = []
    for item in data:
        # Change black (also nearly black shades) to new_color
        if item[0] <= 5 and item[1] <= 5 and item[2] <= 5:
            newData.append(new_color)  # Changing black to new_color
        else:
            newData.append(item)

    image.putdata(newData)
    return image

def invert_image_colors(image_path, output_path, new_color=(30,33,41)):  # Default new color is blue
    with Image.open(image_path) as img:
        if img.format == 'GIF':
            frames = []
            for frame in ImageSequence.Iterator(img):
                inverted_frame = invert_colors(frame)
                color_changed_frame = change_black_to_color(inverted_frame, new_color)
                frames.append(color_changed_frame.convert('P', palette=Image.ADAPTIVE))
            frames[0].save(output_path, save_all=True, append_images=frames[1:], loop=0, duration=img.info.get('duration', 100), optimize=False)
        else:
            inverted_image = invert_colors(img)
            color_changed_image = change_black_to_color(inverted_image, new_color)
            color_changed_image.save(output_path)

def process_folder(folder_path, delete_old):
    if delete_old:
        delete_previous_conversions(folder_path)
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(('colorinverted.png', 'colorinverted.gif')):
            print(f"Skipped {filename}, converted image already exists.")
        elif filename.lower().endswith(('.png', '.gif')):
            base_name, extension = os.path.splitext(filename)
            inverted_filename = f"{base_name}-colorinverted{extension}"
            image_path = os.path.join(folder_path, filename)
            output_path = os.path.join(folder_path, inverted_filename)
            # Only process if delete_old is True or the converted file does not exist
            if delete_old or not os.path.exists(output_path):
                invert_image_colors(image_path, output_path)
                print(f"Processed {filename} -> {inverted_filename}")
            else:
                print(f"Skipped {filename}, converted image already exists.")

# Example usage - replace 'path/to/folder' with the actual folder path
process_folder('./docs/img/', delete_old = False)
