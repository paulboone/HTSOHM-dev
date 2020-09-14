#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

output_path = "comparison_graph1.png"
col = 18

csv_files = ["ms_20_1site", "ms_20_2site", "ms_20_4site", "ms_20_8site", "ms_20_16site", "ms_20_32site", "ms_20_64site", "ms_20_512site"]

ylabel = ""
ymax = None
all_data = []
max_mats = 0

print("loading data...")

for path in csv_files:
    r = np.loadtxt(path + ".csv", delimiter=',', skiprows=1, usecols=col, max_rows=20000)
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
ax.grid(linestyle='-', color='0.8', zorder=0)

ax.set_prop_cycle(color=cm.get_cmap("coolwarm",8)(range(8)))

ax.plot(range(max_mats), np_data)
if ymax:
    ax.axhline(ymax, linestyle="--", lw=3, color="black", label="Max")

# ax.legend(csv_files)
ax.legend(csv_files, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

fig.savefig(output_path, dpi=300)
plt.close(fig)
