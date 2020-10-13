
def choose_parents(num_parents, ids, props, specific_index):
    return [ids[specific_index - 1] for _ in range(num_parents)], [props[specific_index - 1] for _ in range(num_parents)]
