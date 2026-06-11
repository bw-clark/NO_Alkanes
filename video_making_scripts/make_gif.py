
import glob
import re
import imageio.v3 as iio
from PIL import Image, ImageDraw, ImageFont, ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True
fp_in = "*.png"
fp_out = "animated_output.v2.gif"
files = glob.glob(fp_in)
if 'para' in files[0]:
    sorted_files = sorted(files, key=lambda x: float(re.split(r'[a_]', x)[2]))
if 'perp' in files[0]:
    sorted_files = sorted(files, key=lambda x: float(re.split(r'[p_]', x)[2]))
images = []
for filename in sorted_files:
    # Extract angle from filename
    if 'para' in files[0]:
        angle = float(re.split(r'[a_]', filename)[2])
    if 'perp' in files[0]:
        angle = float(re.split(r'[p_]', filename)[2])

    # Read image
    img = Image.fromarray(iio.imread(filename))
    draw = ImageDraw.Draw(img)
    # Use default font
    font = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        size=48
    )
    label = f"{angle:.1f}°"
    # Position in bottom-left corner
    bbox = draw.textbbox((0, 0), label, font=font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    x = 10
    y = img.height - text_h - 10
    # Optional white background for readability
    draw.rectangle(
        [x - 5, y - 5, x + text_w + 5, y + text_h + 5],
        fill="white"
    )
    draw.text((x, y), label, fill="black", font=font)
    images.append(img)
# Convert PIL images back to numpy arrays
images = [image.copy() for image in images]
iio.imwrite(
    fp_out,
    images,
    duration=200,  # milliseconds per frame
    loop=0
)

