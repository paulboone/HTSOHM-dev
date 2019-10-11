#!/usr/bin/env python3

import click
import numpy as np
import matplotlib.pyplot as plt

@click.command()
@click.argument('csv-path', nargs=-1, type=click.Path())
@click.option('-c', '--col', type=int)
@click.option('-o', '--output-path', type=click.Path(), default="comparison_graph.png")
@click.option('--ylabel')
def bin_graph(csv_path, col, output_path, ylabel=""):

    all_data = []
    max_mats = 0

    print("loading data...")
    for path in csv_path:
        r = np.loadtxt(path, delimiter=',', skiprows=1, usecols=col)
        rmax_by_mat = np.maximum.accumulate(r)
        all_data.append(rmax_by_mat)
        max_mats = max(max_mats, rmax_by_mat.size)

    np_data = np.empty((max_mats, len(all_data)))
    for i, dataset in enumerate(all_data):
        np_data[:,i] = dataset

    print("plotting...")
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("# materials")
    ax.set_ylabel(ylabel)

    ax.grid(linestyle='-', color='0.8', zorder=0)
    ax.plot(range(max_mats), np_data)
    ax.legend(csv_path)

    fig.savefig(output_path)
    plt.close(fig)



if __name__ == '__main__':
    bin_graph()
