#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

output_path = "comparison_graph.png"

solid_files = ["ms_20_1site", "ms_20_2site", "ms_20_4site", "ms_20_8site", "ms_20_16site", "ms_20_32site", "ms_20_64site", "ms_20_512site"]
dashed_files = ["random_1site", "random_2site", "random_4site", "random_8site", "random_16site", "random_32site", "random_64site", "random_512site"]
solid_col = 18
dashed_col = 18

# solid_files = ["ms_20_1site", "ms_20_2site", "ms_20_4site", "ms_20_8site"]
# dashed_files = ["ms_20_1-1sites", "ms_20_1-2sites", "ms_20_1-4sites", "ms_20_1-8sites"]
# solid_col = 18
# dashed_col = 13

ylabel = ""
ymax = None


print("loading data...")

def load_data(paths, max_rows, col=18):
    max_mats = 0
    all_data = []
    for path in paths:
        r = np.loadtxt(path + ".csv", delimiter=',', skiprows=1, usecols=col, max_rows=max_rows)
        rmax_by_mat = np.maximum.accumulate(r)
        all_data.append(rmax_by_mat)
        max_mats = max(max_mats, rmax_by_mat.size)

    # copy all data into one array that has normalized dimensions
    np_data = np.empty((max_mats, len(all_data)))
    for i, dataset in enumerate(all_data):
        np_data[:,i] = dataset

    return np_data


solid_data = load_data(solid_files, max_rows=25000, col=solid_col)
dashed_data = load_data(dashed_files, max_rows=25000, col=dashed_col)
max_mats = min(len(solid_data), len(dashed_data))

print("plotting...")
fig = plt.figure(figsize=(5.75,3.75), tight_layout=True)
ax = fig.add_subplot(1, 1, 1)

ax.set_xlabel("# materials")
ax.set_ylabel(ylabel)
ax.grid(linestyle='-', color='0.8', zorder=0)

ax.set_prop_cycle(color=cm.get_cmap("coolwarm",len(solid_files))(range(len(solid_files))))

ax.plot(range(max_mats), solid_data)
ax.plot(range(max_mats), dashed_data, linestyle='--')
if ymax:
    ax.axhline(ymax, linestyle="--", lw=3, color="black", label="Max")

# ax.legend(csv_files)
ax.legend(solid_files + dashed_files, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

fig.savefig(output_path, dpi=300)
plt.close(fig)
