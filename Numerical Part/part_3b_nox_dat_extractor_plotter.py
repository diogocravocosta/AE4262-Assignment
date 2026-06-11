import numpy as np
import matplotlib.pyplot as plt
import os

# =====================================================================
# 1. SETUP
# =====================================================================
# Automatically get the exact folder where this python script is saved
# This ensures it looks in "/Users/hal-9000/.../Numerical Part/"
script_dir = os.path.dirname(os.path.abspath(__file__))

# Map phi to the exact filenames (note the .50 instead of .5)
file_map = {
    0.5: os.path.join(script_dir, 'yi_0.50.dat'),
    0.8: os.path.join(script_dir, 'yi_0.80.dat'),
    1.0: os.path.join(script_dir, 'yi_1.00.dat')
}

mw_NO = 30.006  # Molecular weight of NO in g/mol

phis_plot = []
no_ppm_plot = []

# =====================================================================
# 2. DATA EXTRACTION
# =====================================================================
for phi, filepath in sorted(file_map.items()):
    if not os.path.exists(filepath):
        print(f"WARNING: Could not find '{filepath}'. Skipping phi={phi}...")
        continue
        
    with open(filepath, 'r') as file:
        lines = file.readlines()
        
    # Search for the header line containing column names
    header_idx = -1
    for i, line in enumerate(lines):
        if '[FILE_STRUCTURE_COLUMNS_CONTAINING]' in line:
            header_idx = i + 1
            break
            
    if header_idx != -1:
        # Parse the column names
        header = lines[header_idx].split()
        
        # Get the exact column indices for NO and MeanMass
        idx_NO = header.index('NO')
        idx_MeanMass = header.index('MeanMass')
        
        # Read the very last line of the file (Grid point 300: Flame Exhaust)
        # We skip any completely blank lines at the end of the file
        last_line = []
        for line in reversed(lines):
            if line.strip():
                last_line = line.split()
                break
        
        Y_NO = float(last_line[idx_NO])           # NO Mass Fraction
        M_mix = float(last_line[idx_MeanMass])    # Mixture Mean Molecular Mass
        
        # Convert Mass Fraction to Mole Fraction
        X_NO = Y_NO * (M_mix / mw_NO)
        
        # Convert Mole Fraction to ppm by volume
        NO_ppm = X_NO * 1e6
        
        phis_plot.append(phi)
        no_ppm_plot.append(NO_ppm)
        
        print(f"Successfully processed phi={phi} | Exhaust NO: {NO_ppm:.2f} ppm")

# =====================================================================
# 3. PLOTTING
# =====================================================================
if phis_plot:
    plt.figure(figsize=(8, 5))
    
    plt.plot(phis_plot, no_ppm_plot, '-o', color='#2ca02c', linewidth=2, 
             markersize=8, label='1D Flame Simulation ($CH_4$)')

    # The 'r' before the string stops the \p SyntaxWarning
    #plt.title(r'1D Flame Simulation: NO Emissions vs $\phi$', fontsize=14)
    plt.xlabel(r'Equivalence Ratio, $\phi$', fontsize=12)
    plt.ylabel('NO Emissions (ppm)', fontsize=12)
    
    # Set y-axis to logarithmic scale
    plt.yscale('log')
    plt.ylim([1e-2, 1e4])
    # Ensure the x-axis strictly shows our simulated points
    plt.xticks(phis_plot)
    
    # Turn on major and minor gridlines for logarithmic readability
    plt.grid(True, which='major', linestyle='-', alpha=0.6)
    plt.grid(True, which='minor', linestyle='--', alpha=0.3)
    
    plt.legend(fontsize=12, loc='upper left')
    plt.tight_layout()
    
    # Save and display the figure right next to the script
    plt.savefig(os.path.join(script_dir, 'Part3b_1D_Flame_NO_ppm.pdf'))
    plt.show()
else:
    print("\nNo data was plotted. Please check the filenames in the file_map dictionary.")