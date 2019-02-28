from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, ConnectionPatch
import numpy as np
from scipy.spatial import Delaunay


def delaunay_figure(box_r, convergence_bins, output_path, triang=None, children=[], parents=[],
                    bins=[], new_bins=[], title="", patches=None, prop1range=(0.0,1.0), prop2range=(0.0,1.0),
                    perturbation_methods=None, showgrid=False, showtriangulation=True, bin_saturated=10):

    if not triang:
        triang = Delaunay(box_r)

    hull_point_indices = np.unique(triang.convex_hull.flatten())
    hull_points = np.array([box_r[p] for p in hull_point_indices])

    # plot visualization
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])

    ax.set_xticks(prop1range[1] * np.array(range(0,convergence_bins + 1))/convergence_bins)
    ax.set_yticks(prop2range[1] * np.array(range(0,convergence_bins + 1))/convergence_bins)
    ax.tick_params(labelbottom=False, labelleft=False)
    if showgrid:
        ax.grid(linestyle='-', color='0.5', zorder=0)

    dbinx = prop1range[1] / convergence_bins
    dbiny = prop2range[1] / convergence_bins

    total_materials = bins.sum()
    bin_rects = []
    for b, bcount in np.ndenumerate(bins):
        if bcount > 0.0:
            bin_rects.append(Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny, \
                        facecolor=str(max(1-bcount/bin_saturated, 0.0)), edgecolor='0.8'))
    pc = PatchCollection(bin_rects, match_original=True)
    ax.add_collection(pc)

    new_bin_rects = [Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny) for b in new_bins]
    pc2 = PatchCollection(new_bin_rects, facecolor='#82b7b7')
    ax.add_collection(pc2)



    # plot all points as triangulation
    if showtriangulation:
        ax.triplot(box_r[:,0], box_r[:,1], triang.simplices.copy(), color='#78a7cc', lw=1)

    # plot hull and labels
    ax.plot(hull_points[:,0], hull_points[:,1], color='#78a7cc', marker='o', linestyle='None', zorder=10)
    # for p in hull_points:
    #     ax.annotate(i, (p[0], p[1] - 0.01), zorder=30, ha="center", va="top", size=9)

    # plot children
    ax.plot(children[:,0], children[:,1], color='#81ff6b', marker='o', linestyle='None', zorder=12)

    # plot chosen parents with proportional circles and label
    if parents:
        parent_counter = Counter([tuple(x) for x in parents]) #need tuples because they are hashable
        unique_parents = np.array([[x[0], x[1], num] for x, num in parent_counter.items()])
        ax.scatter(unique_parents[:,0], unique_parents[:,1], s=40*unique_parents[:,2], color='#ffbe6b', marker='o', linestyle='None', zorder=15)
        for p in unique_parents:
            x, y, parent_count = p
            if parent_count > 5:
                ax.annotate(str(int(parent_count)), (x, y), zorder=30, ha="center", va="center", size=9)

    if patches == "donut":
        ax.add_patch(Circle((0.5, 0.5), 0.125, fill=False, linestyle="dashed", linewidth=1, zorder=50))
        ax.add_patch(Circle((0.5, 0.5), 0.375, fill=False, linestyle="dashed", linewidth=2, zorder=50))

    if perturbation_methods:
        p_cmap = {"lattice": "w", "density": "black", "atom_types": "#0FA3B1", "atom_sites": "#F77936", }
        arrows = [ConnectionPatch(parents[i], children[i], "data", arrowstyle="-|>", \
                    shrinkA=5, shrinkB=5, mutation_scale=20, fc=p_cmap[perturbation_methods[i]], \
                    linestyle="--", zorder=40) \
                    for i in range(len(children))]
        for a in arrows:
            ax.add_patch(a)

    ax.set_title(title)
    fig.savefig(output_path)
    plt.close(fig)
