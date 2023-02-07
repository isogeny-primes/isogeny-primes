from sage.all import Integer


def sage_converter(object):
    if isinstance(object, Integer):
        return int(object)

    if isinstance(object, set):
        return sorted(object)

    raise TypeError(f"Object {str(object)} of type {type(object)} is not JSON serializable")
