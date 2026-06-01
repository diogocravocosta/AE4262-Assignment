import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 1. Set up the dimensions (Width, Height) to make it wide and thin
fig, ax = plt.subplots(figsize=(16, 1.5))

# 2. Define the gradient colors: Disagree (Red) to Agree (Green)
# Using hex codes for vibrant, clean colors
colors = ["#e63946", "#0bd430"] 
cmap = LinearSegmentedColormap.from_list("agree_disagree", colors)

# 3. Create the data for the gradient
# A 1D array of 1028 points, stacked so it forms a 2D rectangle
gradient = np.linspace(0, 1, 1028)
gradient = np.vstack((gradient, gradient))

# 4. Plot the gradient
# We use pcolormesh instead of imshow so it saves as true vector geometry, not a rasterized image inside an SVG
ax.pcolormesh(gradient, cmap=cmap, shading='gouraud')

# 5. Remove axes, borders, and ticks for a clean rectangle
ax.set_axis_off()

# 6. Add the text labels inside the rectangle
# 'transform=ax.transAxes' ensures the coordinates (0 to 1) are relative to the rectangle's size
ax.text(0.02, 0.5, 'Disagree', color='white', fontsize=18, fontweight='bold',
        va='center', ha='left', transform=ax.transAxes)

ax.text(0.98, 0.5, 'Agree', color='white', fontsize=18, fontweight='bold',
        va='center', ha='right', transform=ax.transAxes)

# 7. Save as a vector graphic (SVG) and display it
plt.tight_layout()
plt.savefig("agree_disagree_spectrum.svg", format="svg", bbox_inches='tight', pad_inches=0)
plt.show()