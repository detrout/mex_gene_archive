
def cast_to_list(value):
    """Make sure value passed in is a list

    :param value: a value that might be a string or list

    :returns: the list passed in or the value as a single item list
    """
    if isinstance(value, list):
        return value
    elif isinstance(value, str):
        return [value]
