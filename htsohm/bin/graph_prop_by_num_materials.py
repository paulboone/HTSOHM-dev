#!/usr/bin/env python3
import os

import click
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

@click.command()
@click.argument('csv-path', nargs=-1, type=click.Path())
@click.option('-c', '--col', type=int)
@click.option('-o', '--output-path', type=click.Path(), default="comparison_graph.png")
@click.option('--ylabel')
@click.option('--ymax', type=float)
@click.option('--xmax', type=float)
@click.option('--colors', default="Viridis")
def bin_graph(csv_path, col, output_path, ylabel="", ymax=None, xmax=None, colors="Greens"):
    print("loading data...")
    np_data = np.vstack([np.loadtxt(path, delimiter=',', skiprows=1, usecols=col) for path in csv_path])
    legend_labels = [os.path.splitext(os.path.basename(f))[0] for f in csv_path]

    print("plotting...")
    fig = plt.figure(figsize=(5.75,3.75), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("# materials")
    ax.set_ylabel(ylabel)

    # line_colors = cm.get_cmap(colors, len(csv_path))(range(0,len(csv_path)))
    line_colors = cm.get_cmap(colors, len(csv_path) + 3)(range(2,len(csv_path) + 2))

    ax.grid(linestyle='-', color='0.8', zorder=0)
    ax.set_prop_cycle(color=line_colors)

    if xmax:
        ax.set_xlim(0,xmax)

    for i, row in enumerate(np_data):
        # print(row.shape)
        if legend_labels[i] == "reference":
            ax.plot(range(len(row)), row, lw=1, marker="s", markevery=5000, markersize=3)
        else:
            ax.plot(range(len(row)), row, lw=1)

    if ymax:
        ax.axhline(ymax, linestyle="--", lw=2, color="black", label="Max")


    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    fig.savefig(output_path, dpi=300)
    plt.close(fig)



if __name__ == '__main__':
    bin_graph()
