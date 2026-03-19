def get_degree(c: str) -> int:
    """
    Get the degree of the feature based on number of spaces in feature name !?!
    """
    if c == "const":
        return 0
    else:
        return c.count(" ") + 1
