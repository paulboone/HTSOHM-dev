def choose_parents(num_parents, ids, props):
    mats = [(i, m[1]) for i, m in enumerate(props)]
    mats.sort(key=lambda x: x[1])

    # since we are sorted, these are the materials with the highest abs value
    mats = mats[-num_parents:]
    parent_indices = [i for i,m in mats]

    return [ids[i] for i in parent_indices], [props[i] for i in parent_indices]

def choose_specific_parent(num_parents, ids, props, specific_index):
    return [ids[specific_index] for _ in num_parents], [props[specific_index] for _ in num_parents]
