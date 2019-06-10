
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


range_start = np.linspace(0,40, 8, endpoint=False)
dr = range_start[1] - range_start[0]


# pertubation_types = ["all"]
pertubation_types = ["lattice", "lattice_nodens", "atom_types", "atom_sites", "density", "all"]
for t in pertubation_types:
    data = np.load("%s.npy" % t)
    df = pd.DataFrame(data, columns=["ml", "dml"])
    averages = [df[(df.ml < a + dr) & (df.ml >= a)].mean().dml for a in range_start]

    fig = plt.figure(figsize=(5,5), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(range_start + dr/2, averages, zorder=3, color="orange", lw="3")
    ax.scatter(data[:, 0], data[:, 1], zorder=2)
    ax.set_xlabel('parent ML')
    ax.set_ylabel('∆ML')
    # ax.set_xlabel('Parent methane loading [bin units]')
    # ax.set_ylabel('∆ methane loading [bin units]')
    ax.grid(linestyle='-', color='0.7', zorder=0)
    ax.set_xlim(0,40)
    ax.set_ylim(-15,15)
    ax.legend(["Range average"])
    ax.set_title(t)

    fig.savefig("%s.png" % t)
    plt.close(fig)
