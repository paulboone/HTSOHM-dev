#!/usr/bin/env python3
import click
import cmocean
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

prop1range = [-0.01, 1.0]   # max_pair_distance
prop2range = [0.0, 800.0] # ML
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx
fsl = fs = 8

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'weight' : 'bold',
        'size'   : 8}

matplotlib.rc('font', **font)
# rc('text', usetex=True)

@click.command()
@click.argument('csv-path', type=click.File())
def figure4_ml_vs_max_pair_distance(csv_path):
    fig = plt.figure(figsize=(8.5 / 2.54, 8.5 / 2.54))

    cm = matplotlib.cm.get_cmap("viridis")
    # cm = cmocean.cm.thermal
    points = pd.read_csv(csv_path)
    points['ch4_uc'] = points.absolute_volumetric_loading * (num_ch4_a3 * points.a * points.b * points.c)

    ax = fig.subplots(ncols=1)

    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])
    ax.set_xticks(prop1range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    ax.set_yticks(prop2range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))

    ax.tick_params(axis='x', which='major', labelsize=fs)
    ax.tick_params(axis='y', which='major', labelsize=fs)

    ax.grid(which='major', axis='both', linestyle='-', color='0.9', zorder=0)
    # ax.grid(which='minor', axis='both', linestyle='-', color='0.9', zorder=0)

    sc = ax.scatter(points.max_pair_distance, points.absolute_volumetric_loading, zorder=2,
                alpha=0.6, s=points.a, edgecolors=None, linewidths=0, c=points.ch4_uc.round(),
                # c=np.log(points.epsilon_density),
                # norm=matplotlib.colors.LogNorm(vmin=points.epsilon_density.min(), vmax=points.epsilon_density.max()),
                cmap=cm)

    ax.axvline(3**0.5 / 2, 0, 1, lw=1, linestyle="--", color="0.5", label="Site distribution max", zorder=1)

    ax.set_xlabel('Site Distribution', fontsize=fsl)
    ax.set_ylabel('Methane Loading [V/V]', fontsize=fsl)

    # fig.subplots_adjust(wspace=0.05, hspace=0.05)
    output_path = "figure4.png"
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    figure4_ml_vs_max_pair_distance()
