#!/usr/bin/env python3
import os

import click
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# from matplotlib import cm
import pandas as pd

font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'weight' : 'normal',
        'size'   : 8}

rc('font', **font)



@click.command()
def figure_ime_vs_random_efficiency():
    print("loading data...")


    def setup_plot(csv_paths):
        df = pd.DataFrame(data=dict(num=range(1, 500001)))
        ax.set_xlabel("Materials")
        ax.set_yticks([0, 400, 800, 1200, 1600])
        ax.set_xlim(0, 150000)
        ax.set_ylim(0, 1600)
        # ax.axes.yaxis.set_visible(False)
        ax.set_xticks([0, 25000, 50000, 100000, 150000])
        ax.axes.xaxis.set_ticklabels(["0", "25K", "50K", "450K", "500K"])
        ax.grid(linestyle='-', color='0.8', zorder=0, axis="x")
        # ax.axhline(1600, linestyle="--", lw=2, color="black", label=0)
        for path in csv_paths:
            df[path] = pd.read_csv(path, usecols=["unique_bins"])
        return df


    fig = plt.figure(figsize=(3.346, 3.346), tight_layout=True)

    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel("Bins Explored")
    legend_labels = ["IME", "Random"]
    g1 = setup_plot(["reference.csv", "random16_500K.csv"])
    ax.plot(g1.num, g1["reference.csv"], lw=1.5, color="black", zorder=10)
    ax.plot(g1.num[0:50000], g1["random16_500K.csv"][0:50000], lw=1.5, color="orange", zorder=10)
    ax.plot(g1.num[100000:150000], g1["random16_500K.csv"][450000:500000], lw=1.5, color="orange", zorder=10)
    ax.plot([50000, 100000],[541, 730], lw=1.5, color="orange", linestyle="--", zorder=10)
    # ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., facecolor='white', framealpha=1)
    ax.legend(legend_labels, loc='upper right', facecolor='white', framealpha=1)



    ax.axhline(468, lw=1.0, linestyle="--", color="grey", zorder=1)
    ax.axhline(732, lw=1.0, linestyle="--", color="grey", zorder=1)
    ax.axhline(1062, lw=1.0, linestyle="--", color="grey", zorder=1)

    arrow_args = dict(arrowstyle="->")
    ax.annotate("732 bins @ 494956", xy=(494956 - 350000, 732), xycoords="data",
        textcoords="offset points", xytext=(0, 10), horizontalalignment='right', verticalalignment='bottom',
        arrowprops=arrow_args)
    ax.annotate("732 bins @ 4283", xy=(4283, 732), xycoords="data",
        textcoords="offset points", xytext=(13, 10), horizontalalignment='left', verticalalignment='bottom',
        arrowprops=arrow_args)

    ax.annotate("468 bins @ 25000", xy=(25000, 468), xycoords="data",
        textcoords="offset points", xytext=(0, -10), horizontalalignment='left', verticalalignment='top',
        arrowprops=arrow_args)
    ax.annotate("468 bins @ 1786", xy=(1786, 468), xycoords="data",
        textcoords="offset points", xytext=(10, 15), horizontalalignment='left', verticalalignment='bottom',
        arrowprops=arrow_args)

    ax.annotate("1062 bins @ 25000", xy=(25000, 1062), xycoords="data",
        textcoords="offset points", xytext=(0, 10), horizontalalignment='left', verticalalignment='bottom',
        arrowprops=arrow_args)
    fig.savefig("figure_ime_vs_random.png", dpi=1200)
    plt.close(fig)



if __name__ == '__main__':
    figure_ime_vs_random_efficiency()
