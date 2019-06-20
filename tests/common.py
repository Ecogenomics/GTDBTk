def are_files_equal(path_left, path_right, ignore_order):
    """
    Determines if two files are equal based on their contents.

    :param path_left: The left file to compare.
    :param path_right: The right file to compare.
    :param ignore_order: True if the order of lines should be ignored, False otherwise.
    """

    with open(path_left, 'r') as f:
        lines_left = f.readlines()

    with open(path_right, 'r') as f:
        lines_right = f.readlines()

    if len(lines_left) != len(lines_right):
        return False

    if ignore_order:
        return set(lines_left).union(set(lines_right)) == set(lines_left)

    for line_left, line_right in zip(lines_left, lines_right):
        if line_left != line_right:
            return False
    return True
