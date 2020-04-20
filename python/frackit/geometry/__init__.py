from ._geometry import *

# Construct an n-dimensional point
def Point(*args, **kwargs):

    numArgs = len(args) + len(kwargs)
    if numArgs == 0:
        raise Exception("Arguments needed for point construction. For default-" + \
                        "constructible points use the dimension-specific " + \
                        "implementations Point_1, Point_2 or Point_3");

    # single argument construction
    if numArgs == 1:
        # construction from a list
        if isinstance(args[0], list):
            if len(args[0]) == 1: return Point_1(*args, **kwargs)
            elif len(args[0]) == 2: return Point_2(*args, **kwargs)
            else: return Point_3(*args, **kwargs)
        # construction from a single (hopefully) scalar
        else: return Point_1(*args, **kwargs)

    # select point type from the number of arguments
    if numArgs == 2: return Point_2(*args, **kwargs)
    else: return Point_3(*args, **kwargs)
