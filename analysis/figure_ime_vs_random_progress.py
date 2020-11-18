#!/usr/bin/env python3
import click
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, ConnectionPatch
import numpy as np

from htsohm.htsohm_run import calc_bins


def delaunay_figure(ax, convergence_bins, bins=[], prop1range=(0.0,1.0), prop2range=(0.0,1.0), cmname="Greys", alpha=1.0, bin_saturated=20):
    """modified from figures.py"""

    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])

    # ax.axes.xaxis.set_ticklabels([])
    # ax.axes.yaxis.set_ticklabels([])
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    # ax.set_xticks(prop1range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    # ax.set_yticks(prop2range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    # ax.set_xticks(prop1range[1] * np.array(range(0,convergence_bins + 1))/convergence_bins, minor=True)
    # ax.set_yticks(prop2range[1] * np.array(range(0,convergence_bins + 1))/convergence_bins, minor=True)

    # ax.grid(linestyle='-', color='0.8', zorder=0)

    dbinx = prop1range[1] / convergence_bins
    dbiny = prop2range[1] / convergence_bins

    cm = matplotlib.cm.get_cmap(cmname)
    total_materials = bins.sum()
    bin_rects = []
    for b, bcount in np.ndenumerate(bins):
        if bcount > 0.0:
            bin_rects.append(Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny, \
                        facecolor=cm(0.2 + bcount/bin_saturated), edgecolor=None, alpha=alpha))
            # grayscale str(max(0.9-bcount/bin_saturated, 0.0))
    pc = PatchCollection(bin_rects, match_original=True)
    ax.add_collection(pc)

    # new_bin_rects = [Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny) for b in new_bins]
    # pc2 = PatchCollection(new_bin_rects, facecolor='#82b7b7')
    # ax.add_collection(pc2)


@click.command()
@click.argument('csv-paths', nargs=-1)
def figure2_bin_graph(csv_paths=("reference", "random16")):
    prop1range = [0.0, 1.0]   # VF
    prop2range = [0.0, 800.0] # ML
    num_bins = 40

    materials = [50,500, 5000, 25000]

    # plot visualization
    # NOTE: this is set larger than a half column 3.45 inch because of how matplotlib handles whitespace
    # and tight bounding boxes. This seems to be about right to use with labels to the left.
    fig = plt.figure(figsize=(3.8,3.8*2))

    for i, csv_path in enumerate(csv_paths):
        for j, last_material in enumerate(materials):
            bins = np.loadtxt(csv_path + ".csv", delimiter=',', skiprows=1, usecols=(16, 17), max_rows=last_material, dtype=int)
            bin_counts = np.zeros((num_bins, num_bins))
            for bx, by in bins:
                bin_counts[bx,by] += 1
            ax = fig.add_subplot(len(materials),len(csv_paths), len(csv_paths)*j + i + 1)
            for axis in ['top','bottom','left','right']:
              ax.spines[axis].set_linewidth(0.5)
              ax.spines[axis].set_color("gray")
            delaunay_figure(ax, num_bins, bin_counts,prop1range, prop2range)

    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    output_path = "figure2-%s.png" % "-".join(csv_paths)
    fig.savefig(output_path, dpi=1200, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    figure2_bin_graph()
