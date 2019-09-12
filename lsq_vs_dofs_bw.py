import os

import numpy as np
import matplotlib.pyplot as plt


def figure_lsq_vs_dof(data, starting_lsq, output_path, mutation_strength, ymax=0.4):
    # path = "ms20_normal_25.np"
    # filename = os.path.splitext(os.path.basename(path))[0]
    # data = np.loadtxt(path)

    data = data.transpose()

    fig = plt.figure(figsize=(6,6), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    ax.grid(axis="y")
    ax.set_yscale("log")
    ax.set_ylim(starting_lsq*1e-2, starting_lsq*1e+2)
    ax.set_xlabel("DOF multiplier (actual DOF 3x)")
    ax.set_ylabel("LSQ")

    ax.set_title("1000 perturbations, mutation_strength=%3.1f%%" % (mutation_strength * 100))
    ax.boxplot(data, sym='.')

    ax.axhline(starting_lsq, lw=2, linestyle="--", color="black", label="starting lsq = %5.4f" % starting_lsq)
    ax.legend()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
