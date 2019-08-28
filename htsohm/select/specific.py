
def choose_parents(num_parents, box_d, box_range, specific_index):
    return [box_d[specific_index] for _ in range(num_parents)], [box_range[specific_index] for _ in range(num_parents)], None
