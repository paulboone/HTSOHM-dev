
import itertools
from random import uniform
from time import perf_counter

import numpy as np

from lsq_vs_dofs_bw import figure_lsq_vs_dof

def perturb_circular(v, max_change, hlimit = 1.0):
    return (v + uniform(-max_change, max_change)) % hlimit

def deltas(a, hlimit=1.0):
    pairs = itertools.combinations(a, 2)
    # print(list(pairs))
    deltas = []
    for p in pairs:
        p1, p2 = p
        # print(p1, p2)
        d = abs(p2 - p1)
        deltas += [min(d, 1 - d)]
    return deltas

def lsq_error(a, normalization_constant=1):
    return sum([v*v for v in a]) / normalization_constant **2


## WARNING: only good when smallest delta is always within bounds!
def lsq_error_deltas(a, hlimit=1.0, normalization_constant=1):
    pairs = itertools.combinations(a, 2)
    a_np = np.array(list(pairs))
    deltas = a_np[:,1] - a_np[:,0]
    return np.sum(deltas ** 2) / normalization_constant **2

def run_simulations(template_a, max_mult, ms):
    start_time = perf_counter()

    agg_results = np.zeros((max_mult, 1000))
    orig_error = lsq_error_deltas(template_a)

    for dofmult in range(1, max_mult + 1):
        print("------------------------------")
        print("dofmult: %d" % dofmult)
        a = template_a * dofmult
        print(("starting: err %f; mean %f points" % (orig_error, np.mean(a))), a)
        errors = []
        arrays = []
        for _ in range(1000):
            a1 = [perturb_circular(v, ms) for v in a]
            arrays += [a1]
            errors += [lsq_error_deltas(a1, normalization_constant=dofmult)]

        num_better = np.where(np.array(errors) < orig_error)[0]
        print("%d/%d improved: %4.2f" % (len(num_better), len(arrays), len(num_better) / len(arrays)))

        best_error = np.amin(errors)
        best_index = np.where(errors == best_error)[0][0]
        print(("best: err %f; mean %f points" % (best_error, np.mean(arrays[best_index]))), arrays[best_index])
        agg_results[dofmult - 1,:] = errors

    figure_lsq_vs_dof(agg_results, orig_error, "ms%d_%6.5f_%d.png" % (200*ms, orig_error, max_mult), ms*2)
    np.savetxt("ms%d_%6.5f_%d.np" % (200*ms, orig_error, max_mult), agg_results)

    end_time = perf_counter()
    print("Elapsed seconds: %4.2f" % (end_time - start_time))


# lsq_error(deltas([0.29, 0.30, 0.29]), normalization_constant=1)
# lsq_error(deltas([0.29, 0.30, 0.29] * 5), normalization_constant=5)
# lsq_error(deltas([0.25, 0.30, 0.35]))
# lsq_error_deltas([0.20, 0.30, 0.40] * 25, normalization_constant=25)


max_mult = 5

a_bad = [0.25, 0.50, 0.75]
a_normal = [0.20, 0.30, 0.40]
a_better = [0.25, 0.30, 0.35]
a_near_best = [0.29, 0.30, 0.31]
a_best = [0.30, 0.30, 0.30]

templates = [a_bad, a_normal, a_better, a_near_best, a_best]
# templates = [a_near_best, a_best]
mutation_strengths = [0.1, 0.05, 0.025, 0.01]
# mutation_strengths = [0.025]
for t in templates:
    for ms in mutation_strengths:
        run_simulations(t, max_mult, ms)
