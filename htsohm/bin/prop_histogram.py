#!/usr/bin/env python3

import csv

import click
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.argument('prop-column', type=int)
@click.argument('num-bins', type=int)
@click.argument('prop-min', type=int)
@click.argument('prop-max', type=int)
@click.argument('csv-path', type=click.Path())
@click.option('--prop-line', type=float)
def bin_graph(prop_column, num_bins, prop_min, prop_max, csv_path, prop_line=None):

    print("loading materials...")
    csvrows = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    mats = [m[prop_column] for m in csvrows]

    print("plotting...")
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(prop_min, prop_max)
    ax.set_xlabel("Methane Loading (V [STP]/V)")
    ax.set_ylabel("# Materials")

    ax.grid(linestyle='-', color='0.8', zorder=50)

    ax.hist(mats, num_bins, (prop_min, prop_max), rwidth=0.85, alpha=0.9)

    if prop_line:
        ax.axvline(prop_line, 0, 1, lw=2, linestyle="--", color="black", label="Previous high", zorder=100)

    ax.legend()
    ax.set_title("Methane loading of 1000 materials generated from previous high")

    fig.savefig("hist.png")
    plt.close(fig)


if __name__ == '__main__':
    bin_graph()
