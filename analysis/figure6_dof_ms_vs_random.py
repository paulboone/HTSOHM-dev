#!/usr/bin/env python3
import click
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np

from analysis.figure2_algorithm_progress import delaunay_figure

def csv_bin_counts(csv_path, num_bins, max_rows=None):
    bins = np.loadtxt(csv_path  , delimiter=',', skiprows=1, usecols=(16, 17), max_rows=max_rows, dtype=int)
    bin_counts = np.zeros((num_bins, num_bins))
    for bx, by in bins:
        bin_counts[bx,by] += 1
    return bin_counts

def verts_horz_edges(bin_counts):
    """calculate horizontal and vertical lines of random exploration boundary"""
    bins_x, bins_y = bin_counts.shape
    verts = np.zeros((bins_x + 1, bins_y))
    horz = np.zeros((bins_x, bins_y + 1))
    for (x,y), bcount in np.ndenumerate(bin_counts):
        if bcount > 0:
            verts[x, y] += 1
            verts[x + 1, y] += 1
            horz[x, y] += 1
            horz[x, y + 1] += 1
    return verts, horz

def bin_graph_overlay(ax, csv_path_ms, csv_path_random, last_material=None):
    prop1range = [0.0, 1.0]   # VF
    prop2range = [0.0, 800.0] # ML
    num_bins = 40

    bin_counts_ms = csv_bin_counts(csv_path_ms, num_bins, last_material)
    bin_counts_rd = csv_bin_counts(csv_path_random, num_bins, last_material)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)
        ax.spines[axis].set_color("gray")

    delaunay_figure(ax, num_bins, bin_counts_ms, prop1range, prop2range, cmname="Greys", bin_saturated=10)

    # plot boundary of random explored space
    verts, horz = verts_horz_edges(bin_counts_rd)
    dbinx = prop1range[1] / num_bins
    dbiny = prop2range[1] / num_bins
    for (x,y), h in np.ndenumerate(horz):
        if h == 1:
            ax.plot([x*dbinx, (x + 1)*dbinx], [y*dbiny, y*dbiny], lw=1, color="Orange")

    for (x,y), v in np.ndenumerate(verts):
        if v == 1:
            ax.plot([x*dbinx, x*dbinx], [y*dbiny, (y + 1)*dbiny], lw=1, color="Orange")

    # fill random explored space with a transparent green
    bin_rects = []
    for b, bcount in np.ndenumerate(bin_counts_rd):
        if bcount > 0:
            bin_rects.append(Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny, facecolor="Orange", alpha=0.40, fill=True, edgecolor=None))
    pc = PatchCollection(bin_rects, match_original=True)
    ax.add_collection(pc)
    return ax



@click.command()
@click.option('-o', '--output-path', type=click.Path(), default=None)
def figure2_bin_graph(output_path):
    prop1range = [0.0, 1.0]   # VF
    prop2range = [0.0, 800.0] # ML
    num_bins = 40
    csv_paths = [
        ["ms_20_1site.csv", "random_1site.csv"],
        ["ms_20_2site.csv", "random_2site.csv"],
        ["ms_20_4site.csv", "random_4site.csv"],
        ["ms_20_8site.csv", "random_8site.csv"],
        ["ms_20_16site.csv", "random_16site.csv"],
        ["ms_20_32site.csv", "random_32site.csv"],
        ["ms_20_64site.csv", "random_64site.csv"],
        ["ms_20_512site.csv", "random_512site.csv"]
        ]


    # NOTE: this is set larger than a half column 3.45 inch because of how matplotlib handles whitespace
    # and tight bounding boxes. This seems to be about right to use with labels to the left.
    fig = plt.figure(figsize=(4.15 * 2, 4.15))

    for i, (csv_path_ms, csv_path_random) in enumerate(csv_paths):
        ax = fig.add_subplot(2, 4, i + 1)
        bin_graph_overlay(ax, csv_path_ms, csv_path_random, last_material=25000)
        # ax.set_title("512 sites")

    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.savefig(output_path, dpi=1200, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    figure2_bin_graph()
