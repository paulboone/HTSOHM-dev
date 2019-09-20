import os

import numpy as np
import matplotlib.pyplot as plt


def figure_lsq_vs_dof(data, starting_lsq, output_path, mutation_strength, template_max_delta, exp_improvements=None, ymax=0.4):
    # path = "ms20_normal_25.np"
    # filename = os.path.splitext(os.path.basename(path))[0]
    # data = np.loadtxt(path)

    data = data.transpose()
    normalization_factor = 1 / starting_lsq
    data *= normalization_factor

    fig = plt.figure(figsize=(6,6), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    ax.grid(axis="y")
    ax.set_yscale("log")
    ax.tick_params(labelsize=8)

    ax.set_xlabel("DOF multiplier (actual DOF 3x)")
    ax.set_ylabel("LSQnorm / starting LSQnorm")

    ax.set_title("1000 perturbations, mutation_strength=%3.1f%%" % (mutation_strength * 100))
    ax.boxplot(data, sym='.')

    ax.set_ylim(1e-2, 1e+2)
    ax.set_xlim(0, data.shape[1] + 1)

    if starting_lsq != 0.0:
        ax.axhline(starting_lsq * normalization_factor, lw=2, linestyle="--", color="black", label="starting lsq = %5.4f; max ∆ = %4.3f" % (starting_lsq, template_max_delta))

    if exp_improvements:
        for i, improvement in enumerate(exp_improvements):
            ax.text(i + 1, 1.1e-2 + 2e-3 * (i%3), "%6.3f" % abs(100*improvement / starting_lsq), fontsize=6, ha='center')

    ax.legend(loc='upper right')
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
