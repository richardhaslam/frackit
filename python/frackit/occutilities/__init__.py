from ._occutilities import *

###########################################################
# boolean operators for shape wrapper or geometry classes #
###########################################################

# helper function to get the wrapped shape of an object
def getShapeWrapper(object):
    try: return OCCShapeWrapper(_occutilities.getShape(object))
    except: raise TypeError("Given object cannot be wrapped into an occ shape")

# cut operator
def cut(object, tool, eps):
    """
    Boolean operation for cuts.

    Parameters:
    object: The shape to be cut.
    tool: The shape to be cut out of the object.
    eps: Tolerance value to be used.
    """
    return _occutilities.cut(getShapeWrapper(object), getShapeWrapper(tool), eps)

# intersection operator
def intersect(object, tool, eps):
    """
    Boolean operation for intersections.

    Parameters:
    object: The shape whose intersection with the tool is to be computed.
    tool: The tool shape.
    eps: Tolerance value to be used.
    """
    return _occutilities.intersect(getShapeWrapper(object), getShapeWrapper(tool), eps)

# fragmentation operator
def fragment(shapeList, eps):
    """
    Boolean operation for fragmentation.

    Parameters:
    shapeList: List of shapes whose fragments are to be computed.
    eps: Tolerance value to be used.
    """
    return _occutilities.fragment([getShapeWrapper(shape) for shape in shapeList], eps)

# fusion operator
def fuse(shapeList, eps):
    """
    Boolean operation for fusions.

    Parameters:
    shapeList: List of shapes whose fusion is to be computed.
    eps: Tolerance value to be used.
    """
    return _occutilities.fuse([getShapeWrapper(shape) for shape in shapeList], eps)
