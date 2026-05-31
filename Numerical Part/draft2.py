import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================
# 1. TIMELINE DATA
# ==========================================
timeline_data = [
    {"step": "01", "title": "Introduction", "time": "00:00"},
    {"step": "02", "title": "Take a Stand", "time": "07:00"},
    {"step": "03", "title": "Roleplay Simulation", "time": "15:00"},
    {"step": "04", "title": "Reflection and Discussion", "time": "45:00"}
]

# Single color matching the hue of your presentation slides
theme_blue = '#004F91'

# ==========================================
# 2. SETUP THE FIGURE (Transparent Background)
# ==========================================
fig, ax = plt.subplots(figsize=(14, 5))
ax.set_axis_off() 

# Make both the figure and axes backgrounds entirely transparent
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

# Metrics for drawing
num_steps = len(timeline_data)
box_w = 1.0 / num_steps  
box_h = 0.15             
point_w = box_w * 0.2    
start_y = 0.5            

# ==========================================
# 3. DRAW THE TIMELINE
# ==========================================
for i, data in enumerate(timeline_data):
    x = i * box_w
    
    if i == 0:
        verts = [
            (x, start_y), 
            (x + box_w - point_w, start_y), 
            (x + box_w, start_y + box_h/2), 
            (x + box_w - point_w, start_y + box_h), 
            (x, start_y + box_h)
        ]
        # Visual center for the first item (accounts for the flat left edge)
        center_x = x + (box_w - point_w) / 2
    else:
        verts = [
            (x, start_y), 
            (x + box_w - point_w, start_y), 
            (x + box_w, start_y + box_h/2), 
            (x + box_w - point_w, start_y + box_h), 
            (x, start_y + box_h),
            (x + point_w, start_y + box_h/2) 
        ]
        # True visual center for indented items (the indent and point perfectly balance out)
        center_x = x + (box_w / 2)
        
    poly = patches.Polygon(verts, facecolor=theme_blue, edgecolor='white', linewidth=2, zorder=2)
    ax.add_patch(poly)
    
    # Text inside the chevrons
    ax.text(center_x, start_y + box_h/2 + 0.015, data["step"], color='white', 
            fontsize=12, fontweight='bold', ha='center', va='center', zorder=3)
    ax.text(center_x, start_y + box_h/2 - 0.025, data["title"], color='white', 
            fontsize=10, fontweight='bold', ha='center', va='center', zorder=3)

    # ==========================================
    # 4. DRAW THE PINS AND TIMES
    # ==========================================
    pin_length = 0.15
    if i % 2 == 0:
        # Point UP
        line_y_start = start_y + box_h
        line_y_end = line_y_start + pin_length
        
        # Pin line and dots
        ax.plot([center_x, center_x], [line_y_start, line_y_end], color='gray', linewidth=1.5, zorder=1)
        ax.scatter([center_x], [line_y_end], color=theme_blue, s=120, edgecolor='gray', linewidth=1.5, zorder=3)
        ax.scatter([center_x], [line_y_end], color='white', s=30, zorder=4) 
        
        # Time Text (Now White)
        ax.text(center_x, line_y_end + 0.04, data["time"], color='white', 
                fontsize=11, fontweight='bold', ha='center', va='bottom')
    else:
        # Point DOWN
        line_y_start = start_y
        line_y_end = line_y_start - pin_length
        
        # Pin line and dots
        ax.plot([center_x, center_x], [line_y_start, line_y_end], color='gray', linewidth=1.5, zorder=1)
        ax.scatter([center_x], [line_y_end], color=theme_blue, s=120, edgecolor='gray', linewidth=1.5, zorder=3)
        ax.scatter([center_x], [line_y_end], color='white', s=30, zorder=4) 
        
        # Time Text (Now White)
        ax.text(center_x, line_y_end - 0.04, data["time"], color='white', 
                fontsize=11, fontweight='bold', ha='center', va='top')

# ==========================================
# 5. RENDER AND SAVE
# ==========================================
plt.tight_layout()

# transparent=True completely removes the white background
plt.savefig("activity_timeline_transparent_white_text.svg", format="svg", bbox_inches='tight', pad_inches=0.1, transparent=True)
plt.show()