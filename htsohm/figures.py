from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle
import numpy as np
from scipy.spatial import Delaunay


def delaunay_figure(box_r, convergence_bins, output_path, triang=None, children=[], parents=[],
                    bins=[], new_bins=[], title="", patches=None, prop1range=(0.0,1.0), prop2range=(0.0,1.0)):

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
    ax.grid(linestyle='-', color='0.5', zorder=0)

    dbinx = prop1range[1] / convergence_bins
    dbiny = prop2range[1] / convergence_bins
    bin_rects = [Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny) for b in bins]
    pc = PatchCollection(bin_rects, facecolor='0.75')
    ax.add_collection(pc)

    new_bin_rects = [Rectangle((b[0] * dbinx, b[1] * dbiny), dbinx, dbiny) for b in new_bins]
    pc2 = PatchCollection(new_bin_rects, facecolor='#82b7b7')
    ax.add_collection(pc2)

    # plot all points as triangulation
    ax.triplot(box_r[:,0], box_r[:,1], triang.simplices.copy(), color='#78a7cc', lw=1)

    # plot hull and labels
    ax.plot(hull_points[:,0], hull_points[:,1], color='#78a7cc', marker='o', linestyle='None', zorder=10)
    # for p in hull_points:
    #     ax.annotate(i, (p[0], p[1] - 0.01), zorder=30, ha="center", va="top", size=9)

    # plot children
    ax.plot(children[:,0], children[:,1], color='#81ff6b', marker='o', linestyle='None', zorder=12)

    # plot chosen parents with proportional circles and label
    if len(parents) > 0:
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

    ax.set_title(title)
    fig.savefig(output_path)
    plt.close(fig)
