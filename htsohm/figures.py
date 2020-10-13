from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, ConnectionPatch
import numpy as np
from scipy.spatial import Delaunay


def delaunay_figure(props, convergence_bins, output_path, triang=None, children=[], parents=[],
                    bins=[], new_bins=[], title="", patches=None, prop1range=(0.0,1.0),
                    prop2range=(0.0,1.0),perturbation_methods=None, show_grid=True,
                    show_triangulation=True, show_hull=True, bin_saturated=10, bin_scores=None,
                    addl_data_set=None):


    # plot visualization
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])
    ax.set_xlabel("Void Fraction")
    ax.set_ylabel("Henry's Coefficient [mol / (m3 pa)]")

    ax.set_xticks((prop1range[1] - prop1range[0]) * np.array([0.0, 0.25, 0.5, 0.75, 1.0]) + prop1range[0])
    ax.set_yticks((prop2range[1] - prop2range[0]) * np.array([0.0, 0.25, 0.5, 0.75, 1.0]) + prop2range[0])
    ax.set_xticks((prop1range[1] - prop1range[0]) * np.array(range(0,convergence_bins + 1))/convergence_bins + prop1range[0], minor=True)
    ax.set_yticks((prop2range[1] - prop2range[0]) * np.array(range(0,convergence_bins + 1))/convergence_bins + prop2range[0], minor=True)

    if show_grid:
        ax.grid(linestyle='-', color='0.8', zorder=0)

    dbinx = (prop1range[1] - prop1range[0]) / convergence_bins
    dbiny = (prop2range[1] - prop2range[0]) / convergence_bins

    total_materials = bins.sum()
    bin_rects = []
    for b, bcount in np.ndenumerate(bins):
        if bcount > 0.0:
            bin_rects.append(Rectangle((b[0] * dbinx + prop1range[0],
                        prop2range[0] + b[1] * dbiny), dbinx, dbiny,
                        facecolor=str(max(1-bcount/bin_saturated, 0.0)), edgecolor='0.8'))
    pc = PatchCollection(bin_rects, match_original=True)
    ax.add_collection(pc)

    new_bin_rects = [Rectangle((b[0] * dbinx + prop1range[0],
        b[1] * dbiny + prop2range[0]), dbinx, dbiny) for b in new_bins]
    pc2 = PatchCollection(new_bin_rects, facecolor='#82b7b7')
    ax.add_collection(pc2)

    if bin_scores is not None:
        for (i,j), score in bin_scores:
            ax.annotate(str(int(score)), ((i + 0.5) * dbinx, (j+0.5) * dbiny), zorder=80, ha="center", va="center", size=9)

    if not triang and (show_hull or show_triangulation):
        triang = Delaunay(props)

    # plot all points as triangulation
    if show_triangulation:
        ax.triplot(props[:,0], props[:,1], triang.simplices.copy(), color='#78a7cc', lw=1)

    # plot hull and labels
    if show_hull:
        hull_point_indices = np.unique(triang.convex_hull.flatten())
        hull_points = np.array([props[p] for p in hull_point_indices])
        ax.plot(hull_points[:,0], hull_points[:,1], color='#78a7cc', marker='o', linestyle='None', zorder=10)

    # plot children
    if len(children) > 0:
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

    # plot additional data set
    if addl_data_set is not None:
        ax.plot(addl_data_set[:,0], addl_data_set[:,1], color='#FF8C00', marker=',', linestyle='None', zorder=20)


    if patches == "donut":
        ax.add_patch(Circle((0.5, 0.5), 0.125, fill=False, linestyle="dashed", linewidth=1, zorder=50))
        ax.add_patch(Circle((0.5, 0.5), 0.375, fill=False, linestyle="dashed", linewidth=2, zorder=50))

    if perturbation_methods and len(parents) > 0:
        p_cmap = {"all": "black", "lattice": "w", "lattice_nodens": "#848C8E", "density": "#435058", "atom_types": "#0FA3B1", "atom_sites": "#F77936", }
        arrows = [ConnectionPatch(parents[i], children[i], "data", arrowstyle="-|>", \
                    shrinkA=5, shrinkB=5, mutation_scale=20, fc=p_cmap[perturbation_methods[i]], \
                    linestyle="--", zorder=40) \
                    for i in range(len(children))]
        for a in arrows:
            ax.add_patch(a)

    ax.set_title(title)
    fig.savefig(output_path)
    plt.close(fig)
