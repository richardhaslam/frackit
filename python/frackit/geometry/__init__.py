from ._geometry import *

def raiseGeometryConstructorException(type, issue=""):
    if issue == "notImplemented":
        raise NotImplementedError("Requested specialization of " + type + " class is not yet implemented")
    if issue == "numArgs":
        raise RuntimeError("Wrong number of arguments provided for construction of '" + type + "'")
    else:
        raise RuntimeError("Could not construct '" + type + "' from provided argument(s)")


############################################
# Argument-dependent n-d point construction"
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
                        "constructible vectors use the dimension-specific " + \
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
            except: raiseGeometryConstructorException("vector")
        return makeVector(dim)

    # (maybe) construct from two points
    if numArgs == 2:
        try: dim = args[0].worldDimension
        except: pass
        else: return makeVector(dim)

    # last option: construction from the raw coordinates
    return makeVector(numArgs)


############################################
# Argument-dependent n-d direction construction
def Direction(*args, **kwargs):
    numArgs = len(args) + len(kwargs)
    if numArgs != 1: raiseGeometryConstructorException("direction", "numArgs")
    try: dim = args[0].worldDimension
    except: raiseGeometryConstructorException("direction")

    if dim == 1: return Direction_1(*args, **kwargs)
    elif dim == 2: return Direction_2(*args, **kwargs)
    else: return Direction_3(*args, **kwargs)


############################################
# Argument-dependent n-d circle construction
def Circle(*args, **kwargs):
    numArgs = len(args) + len(kwargs)
    if numArgs != 3: raiseGeometryConstructorException("circle", "numArgs")
    try: dim = args[0].worldDimension
    except: raiseGeometryConstructorException("circle")

    if dim == 1: raiseGeometryConstructorException("circle", "notImplemented") # todo
    elif dim == 2: raiseGeometryConstructorException("circle", "notImplemented") # todo
    else: return Circle_3(*args, **kwargs)

############################################
# Compute the distance between two geometries
def computeDistance(geo1, geo2):

    """
    Compute the minimum euclidian distance between geometries.
    If the second argument (geo2) is a list of geometries, the minimum
    distance of geo1 to the geometries of geo2 is returned.
    """

    def doComputation(geo):
        from frackit.occutilities import getShape, OCCShapeWrapper
        try: return _geometry.computeDistance(geo1, geo)
        except: return _geometry.computeDistance(OCCShapeWrapper(getShape(geo1)),
                                                 OCCShapeWrapper(getShape(geo)))

    if not isinstance(geo2, list): return doComputation(geo2)
    else: return min([doComputation(g) for g in geo2])
