import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import time
from math import sqrt

in_file = "C:/Users/Aakanksha Choudhary/OneDrive/Desktop/cosmx/centroids1.csv"
out_file= "C:/Users/Aakanksha Choudhary/OneDrive/Desktop/cosmx/circles/INC2842_C3/my_regions"
sample = "INC2842_C3" 
x_min = 0
x_max = 10000
y_min = 0
y_max = 10000

# read file
df = pd.read_csv(in_file, sep='\t')


# time stamp
time_stamp = time.time()

# subset on sample
df = df[df["SAMPLE"] == sample]

# subset on FOV
df = df[df["X"] >= x_min]
df = df[df["X"] <= x_max]
df = df[df["Y"] >= y_min]
df = df[df["Y"] <= y_max]

# filter columns
df = df.drop('SAMPLE', axis=1)

# get x, y and cell
grid_x = np.array(df["X"])
grid_y = np.array(df["Y"])
cell = np.array(df["CELL_TYPE"])

# colour map
cmap = matplotlib.colors.ListedColormap(["#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"], name='from_list', N=None)

# plot
fig, ax = plt.subplots()
pts = ax.scatter(grid_x, grid_y, c = cell, s=2, cmap=cmap, )

ax.set_facecolor('grey')

def conc_circ(center_x, center_y, x_min, x_max, y_min, y_max, spacing=0.1, **kwargs):
    max_radius = sqrt(max(abs(center_x - x_min), abs(center_x - x_max))**2 +
                  max(abs(center_y - y_min), abs(center_y - y_max))**2)

    radius = spacing
    while radius <= max_radius:
        circle = plt.Circle(
            xy=(center_x, center_y),
            radius=radius,
            fill=False,
            linestyle='--',
            edgecolor='black',
            linewidth=1.5,
            **kwargs
        )
        ax.add_patch(circle)
        radius += spacing

center_x = 2.7
center_y = 2.8

x_min = df['X'].min()
x_max = df['X'].max()
y_min = df['Y'].min()
y_max = df['Y'].max()


conc_circ(center_x=center_x, center_y=center_y, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, spacing=0.1)
plt.show()

ring_cells_dict = {}
spacing = 0.1  # must match conc_circ spacing

# Calculate max_radius again (same formula as in conc_circ)
max_radius = sqrt(max(abs(center_x - x_min), abs(center_x - x_max))**2 +
                  max(abs(center_y - y_min), abs(center_y - y_max))**2)

# Generate radii list
radii = []
radius = spacing
while radius <= max_radius:
    radii.append(radius)
    radius += spacing

print(f"Radii: {radii}")

# Initialize dictionary with ring keys
for i in range(len(radii)):
    ring_key = f'ring_{i+1}'
    ring_cells_dict[ring_key] = []

# Assign cells to rings based on distance
for i, row in df.iterrows():
    X = row['X']
    Y = row['Y']
    cell_id = i

    dist = np.sqrt((X - center_x)**2 + (Y - center_y)**2)

    for j in range(len(radii)):
        inner_radius = 0 if j == 0 else radii[j - 1]
        outer_radius = radii[j]
        if inner_radius < dist <= outer_radius:
            ring_key = f'ring_{j+1}'
            ring_cells_dict[ring_key].append(str(cell_id))
            break

# Output cells per ring
for ring_key, cells in ring_cells_dict.items():
    print(f"{ring_key} contains {len(cells)} cells")
    ring_out_file = (
        f'{out_file}_{sample}_{ring_key}_'
        f'{x_min}_{x_max}_{y_min}_{y_max}_{int(time_stamp)}.txt'
    )
    with open(ring_out_file, 'w') as f:
        for cell in cells:
            f.write(cell + '\n')
