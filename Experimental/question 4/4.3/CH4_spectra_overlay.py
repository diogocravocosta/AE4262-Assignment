from PIL import Image, ImageDraw, ImageFont
import numpy as np

paths = {
    'phi1':  'C:\\Users\\franc\\Documents\\GitHub\\AE4262-Assignment\\Experimental\\question 4\\Spectrometer\\CH4_PHI1.png',
    'phi08': 'C:\\Users\\franc\\Documents\\GitHub\\AE4262-Assignment\\Experimental\\question 4\\Spectrometer\\CH4_PHI0.8.png',
    'phi12': 'C:\\Users\\franc\\Documents\\GitHub\\AE4262-Assignment\\Experimental\\question 4\\Spectrometer\\CH4_PHI1.2.png',
}

tints = {
    'phi12': (0,   80,  200),  # blue  - opaque, back
    'phi1':  (210,  30,   30), # red   - transparent, middle
    'phi08': (0,    0,    0),  # black - opaque, front
}

# phi12 back (opaque), phi1 middle (transparent), phi08 front (opaque)
draw_order = ['phi12', 'phi1', 'phi08']

arrays = {}
sizes = []
for k, p in paths.items():
    img = Image.open(p).convert('L')
    sizes.append(img.size)
    arrays[k] = np.array(img, dtype=np.float32)

W, H = sizes[0]  # use actual image size
result = np.ones((H, W, 3), dtype=np.float32) * 255.0

for k in draw_order:
    arr = arrays[k]
    r, g, b = tints[k]
    signal = 1.0 - arr / 255.0
    alpha = 0.5 if k == 'phi1' else 1.0  # only phi1 is transparent
    effective = signal * alpha
    for c_idx, c_val in enumerate([r, g, b]):
        result[:, :, c_idx] = result[:, :, c_idx] * (1 - effective) + c_val * effective

result = np.clip(result, 0, 255).astype(np.uint8)
out = Image.fromarray(result, 'RGB')

# Resize to 2048 wide keeping aspect ratio
scale = 2048 / W
out_small = out.resize((2048, int(H * scale)), Image.LANCZOS)

# --- Legend ---
draw = ImageDraw.Draw(out_small)
legend_items = [
    ('phi = 0.8', (0,   0,   0)),
    ('phi = 1.0', (210, 30,  30)),
    ('phi = 1.2', (0,   80,  200)),
]

lx, ly = int(out_small.width * 0.75), 80
pad = 18
swatch = 38
line_h = 46
box_w = 320
box_h = line_h * len(legend_items) + 10

draw.rectangle([lx-pad, ly-pad, lx+box_w+pad, ly+box_h+pad], fill=(255,255,255), outline=(160,160,160), width=3)

try:
    font = ImageFont.truetype("C:\\Windows\\Fonts\\arialbd.ttf", 34)
except:
    font = ImageFont.load_default()

for i, (label, fill_color) in enumerate(legend_items):
    y = ly + i * line_h
    border_col = (80, 80, 80) if fill_color == (0, 0, 0) else fill_color
    draw.rectangle([lx, y+4, lx+swatch, y+swatch], fill=fill_color, outline=border_col, width=2)
    draw.text((lx+swatch+14, y), label, fill=(30,30,30), font=font)

out_small.save('CH4_spectra_overlay.png', optimize=True)
print("Done - saved as CH4_spectra_overlay.png")