import itertools
from random import uniform

import numpy as np

def perturb_circular(v, max_change, hlimit = 1.0):
    return (v + uniform(-max_change, max_change)) % hlimit

def deltas(a, hlimit=1.0):
    pairs = itertools.combinations(a, 2)
    # print(list(pairs))
    deltas = []
    for p in pairs:
        p1, p2 = p
        print(p1, p2)
        d = abs(p2 - p1)
        deltas += [min(d, 1 - d)]
    return deltas

def lsq_error(a):
    return sum([v*v for v in a])

a = [0.2, 0.4, 0.3]
orig_error = lsq_error(a)
ms = 0.1 # mutation strength

errors = []
arrays = []
for _ in range(1000):
    a1 = [perturb_circular(v, ms) for v in a]
    arrays += [a1]
    errors += [lsq_error(a1)]

n = np.array([list(range(len(errors))), errors])
num_better = np.where(np.array(errors) < orig_error)[0]
print("%d/%d improved: %4.2f" % (len(num_better), len(arrays), len(num_better) / len(a)))

best_error = np.amin(errors)
best_index = np.where(errors == best_error)[0][0]
print(("best: %f " % best_error), arrays[best_index])
