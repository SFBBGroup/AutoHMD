import numpy as np
import matplotlib.pyplot as plt

# List of data files
files = ['rep0/analysis/iRMSD_r0.dat', 'rep1/analysis/iRMSD_r0.dat', 'rep2/analysis/iRMSD_r0.dat']
colors = ['darkorange', 'darkviolet', 'steelblue']
labels = ['Replica 0', 'Replica 1', 'Replica 2']

# Initialize lists to store RMSD
irmsd_list = []

# Read the data
for file in files:
    data = np.loadtxt(file, skiprows=1)
    irmsd_list.append(data[:, 1])

# Simulation parameters
total_frames = 5601
total_time_ns = 70.0

# Calculate time in nanoseconds for each frame
frame = np.arange(total_frames)
time_ns = (frame / total_frames) * total_time_ns

# Configure legend font
plt.rcParams['legend.fontsize'] = 13

# Create plot
plt.figure(figsize=(10, 6))

# Plot the RMSDs
for i, irmsd in enumerate(irmsd_list):
    plt.plot(time_ns[:len(irmsd)], irmsd, label=labels[i], color=colors[i])

# Add vertical and horizontal lines
plt.axvline(x=30, color='black', linestyle='--')
plt.axvline(x=42.5, color='black', linestyle='--')
plt.axvline(x=55, color='black', linestyle='--')
plt.axhline(y=5, color='black', linestyle='--')

# Add text labels
positions = [15, 36, 48, 62]
temperature_labels = ['310 K', '330 K', '360 K', '390 K']
for pos, temp_label in zip(positions, temperature_labels):
    plt.text(pos, 11.5, temp_label, color='black', verticalalignment='top', 
             horizontalalignment='center', fontsize=22)

# Configure x and y axes
plt.xlabel('Time (ns)', fontsize=24, fontweight='bold')
plt.ylabel('iRMSD (Å)', fontsize=24, fontweight='bold')

# Configure axis intervals
interval_ns = 5.0
tick_markers_ns = np.arange(0, total_time_ns + interval_ns, interval_ns)
plt.xticks(tick_markers_ns, fontsize=22, rotation=45)
plt.yticks(fontsize=22)
plt.ylim(0, 12)
plt.xlim(0, total_time_ns)

# Add legend outside the plot, at the top
plt.legend(loc='upper right', bbox_to_anchor=(0.82, 1.20), ncol=len(labels), fontsize=16)
plt.tight_layout()

plt.savefig(
        fname='iRMSD_replicas.png',  # File name must be a valid string
        transparent=True, 
        dpi=400,  # dpi of 400 for high resolution
        format='png',  # Specifies the image format
        metadata=None, 
        bbox_inches=None, 
        pad_inches=0.5,
        facecolor='auto', 
        edgecolor='auto', 
        backend=None,
       )
plt.grid(False)
plt.show()
