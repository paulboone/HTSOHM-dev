#!/usr/bin/env python3
import click
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

@click.command()
@click.argument('csv-path', type=click.File())
def figure3_ml_vs_vf(csv_path):

    prop1range = [0.0, 1.0]   # VF
    prop2range = [0.0, 800.0] # ML
    num_bins = 40

    num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx

    fsl = fs = 8
    fig = plt.figure(figsize=(4,4))

    cm = matplotlib.cm.get_cmap("viridis")
    points = pd.read_csv(csv_path)
    points['ch4_uc'] = points.absolute_volumetric_loading * (num_ch4_a3 * points.a * points.b * points.c)

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # rc('text', usetex=True)

    ax = fig.subplots(ncols=1)

    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])
    ax.set_xticks(prop1range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    ax.set_yticks(prop2range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    ax.set_xticks(prop1range[1] * np.array(range(0,num_bins + 1))/num_bins, minor=True)
    ax.set_yticks(prop2range[1] * np.array(range(0,num_bins + 1))/num_bins, minor=True)

    ax.tick_params(axis='x', which='major', labelsize=fs)
    ax.tick_params(axis='y', which='major', labelsize=fs)

    ax.set_ylabel("Fraction of applied heat flux", fontsize=fsl)


    ax.grid(which='major', axis='both', linestyle='-', color='0.9', zorder=0)
    # ax.grid(which='minor', axis='both', linestyle='-', color='0.9', zorder=0)

    sc = ax.scatter(points.void_fraction_geo, points.absolute_volumetric_loading, zorder=2,
                alpha=0.6, s=points.a, edgecolors=None, linewidths=0, c=points.ch4_uc.round(),
                cmap=cm)


    ax.set_xlabel('Void Fraction', fontsize=fsl)
    ax.set_ylabel('Methane Loading [V/V]', fontsize=fsl)

    # fig.subplots_adjust(wspace=0.05, hspace=0.05)
    output_path = "figure3.png"
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    figure3_ml_vs_vf()
