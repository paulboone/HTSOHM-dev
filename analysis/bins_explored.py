#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

output_path = "comparison_graph1.png"
col = 18
csv_files = ["reference", "random4", "random8", "random16", "atoms_1", "atoms_2", "atoms_8",
                "ms_5", "ms_10", "ms_40", "types_1", "uff_25", "uff_50"]

# csv_files = ["reference", "atoms_1", "atoms_2", "atoms_8"]
# csv_files = ["reference", "ms_5", "ms_10", "ms_40"]
# csv_files = ["reference", "types_1", "types_4"]


ylabel = ""
ymax = None
all_data = []
max_mats = 0

print("loading data...")

for path in csv_files:
    r = np.loadtxt(path + ".csv", delimiter=',', skiprows=1, usecols=col)
    rmax_by_mat = np.maximum.accumulate(r)
    all_data.append(rmax_by_mat)
    max_mats = max(max_mats, rmax_by_mat.size)

# copy all data into one array that has normalized dimensions
np_data = np.empty((max_mats, len(all_data)))
for i, dataset in enumerate(all_data):
    np_data[:,i] = dataset

print("plotting...")
fig = plt.figure(figsize=(5.75,3.75), tight_layout=True)
ax = fig.add_subplot(1, 1, 1)

ax.set_xlabel("# materials")
ax.set_ylabel(ylabel)


print(cm.get_cmap("Greens", 3)(range(3)))
# colors = [[0, 0, 0, 1]]

colors = [[0, 0, 0, 1]] + \
            cm.get_cmap("Greys", 6)(range(2,5)).tolist() + \
            cm.get_cmap("Greens", 6)(range(2,5)).tolist() + \
            cm.get_cmap("Blues", 6)(range(2,5)).tolist() + \
            cm.get_cmap("Oranges", 6)(range(2,5)).tolist()
            # cm.get_cmap("Reds", 6)(range(2,4)).tolist()
print(colors)
ax.grid(linestyle='-', color='0.8', zorder=0)
# ax.plot(range(max_mats), np_data)
ax.set_prop_cycle(color=colors)
# ax.set_prop_cycle(color=["red", "green", "blue"])
ax.plot(range(max_mats), np_data)
if ymax:
    ax.axhline(ymax, linestyle="--", lw=3, color="black", label="Max")

# ax.legend(csv_files)
ax.legend(csv_files, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

fig.savefig(output_path, dpi=300)
plt.close(fig)
