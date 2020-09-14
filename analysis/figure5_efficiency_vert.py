#!/usr/bin/env python3
import os

import click
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# from matplotlib import cm
import pandas as pd

# fsl = fs = 8
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'weight' : 'normal',
        'size'   : 8}

rc('font', **font)



@click.command()
def figure5_efficiency():
    print("loading data...")


    def setup_plot(csv_paths):
        df = pd.DataFrame(data=dict(num=range(1, 25051)))
        ax.set_xlabel("Materials")
        ax.set_ylabel("Bins Explored")
        ax.set_yticks([0, 400, 800, 1200, 1600])
        ax.set_xlim(0, 25000)
        ax.set_ylim(0, 1600)
        ax.axes.xaxis.set_ticklabels(["0", "5K", "10K", "15K", "20K", "25K"])
        ax.grid(linestyle='-', color='0.8', zorder=0)
        # ax.axhline(1600, linestyle="--", lw=2, color="black", label=0)
        for path in csv_paths:
            df[path] = pd.read_csv(path, usecols=["unique_bins"])

        return df


    # fig = plt.figure(figsize=(6.69,2.0), tight_layout=True)
    fig = plt.figure(figsize=(3.75,8.25), tight_layout=True)

    ax = fig.add_subplot(4, 1, 1)
    legend_labels = ["1", "1-2", "1-4"]
    g1 = setup_plot(["types_1.csv", "reference.csv", "types_4.csv"])
    ax.plot(g1.num, g1["types_1.csv"],      lw=1 )
    ax.plot(g1.num, g1["reference.csv"], lw=1, color="black")
    ax.plot(g1.num, g1["types_4.csv"],     lw=1 )
    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., facecolor='white', framealpha=1)
    ax.set_title("Atom Types")

    ax = fig.add_subplot(4, 1, 2)
    legend_labels = ["1", "1-2", "1-4", "1-8"]
    g1 = setup_plot(["atoms_1.csv", "atoms_2.csv", "reference.csv", "atoms_8.csv"])
    ax.plot(g1.num, g1["atoms_1.csv"],      lw=1 )
    ax.plot(g1.num, g1["atoms_2.csv"],     lw=1 )
    ax.plot(g1.num, g1["reference.csv"], lw=1, color="black")
    ax.plot(g1.num, g1["atoms_8.csv"],     lw=1 )
    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., facecolor='white', framealpha=1)
    ax.set_title("Atom Sites")

    ax = fig.add_subplot(4, 1, 3)
    legend_labels = ["5%", "10%", "20%", "40%", "20% alt", "random"]
    g1 = setup_plot(["ms_5.csv", "ms_10.csv", "reference.csv", "ms_40.csv", "ms_20_questionable_perturbation.csv", "random8.csv"])
    ax.plot(g1.num, g1["ms_5.csv"],      lw=1 )
    ax.plot(g1.num, g1["ms_10.csv"],     lw=1 )
    ax.plot(g1.num, g1["reference.csv"], lw=1, color="black")
    ax.plot(g1.num, g1["ms_40.csv"],     lw=1 )
    ax.plot(g1.num, g1["ms_20_questionable_perturbation.csv"], lw=1, color="pink")
    ax.plot(g1.num, g1["random8.csv"],    lw=1, color="grey" )
    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., facecolor='white', framealpha=1)
    ax.set_title("Mutation Strength")

    ax = fig.add_subplot(4, 1, 4)
    legend_labels = ["UFF", "UFF±25%", "UFF±50%"]
    g1 = setup_plot(["reference.csv", "uff_25.csv", "uff_50.csv"])
    ax.plot(g1.num, g1["reference.csv"], lw=1, color="black")
    ax.plot(g1.num, g1["uff_25.csv"],     lw=1 )
    ax.plot(g1.num, g1["uff_50.csv"],     lw=1 )
    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., facecolor='white', framealpha=1)
    ax.set_title("Sigma / Epsilon")

    fig.savefig("figure5-vert.png", dpi=1200)
    plt.close(fig)



if __name__ == '__main__':
    figure5_efficiency()
