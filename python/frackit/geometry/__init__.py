from ._geometry import *

############################################
# Argument-dependent n-d point construction
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


############################################
# Argument-dependent n-d vector construction
def Vector(*args, **kwargs):

    numArgs = len(args) + len(kwargs)
    if numArgs == 0:
        raise Exception("Arguments needed for vector construction. For default-" + \
                        "constructible vector use the dimension-specific " + \
                        "implementations Vector_1, Vector_2 or Vector_3");

    def makeVector(dim):
        if dim == 1: return Vector_1(*args, **kwargs)
        elif dim == 2: return Vector_2(*args, **kwargs)
        else: return Vector_3(*args, **kwargs)

    # single argument construction (either list or direction)
    if numArgs == 1:
        if isinstance(args[0], list):
            dim = len(args[0])
        else:
            try: dim = args[0].worldDimension
            except: raise Exception("Invalid argument provided for Vector construction")
        return makeVector(dim)

    # (maybe) construct from two points
    if numArgs == 2:
        try: dim = args[0].worldDimension
        except: pass
        else: return makeVector(dim)

    # last option: construction from the raw coordinates
    return makeVector(numArgs)
