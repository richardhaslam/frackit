from ._occutilities import *

# functions of this module provide operations on wrappers around
# OpenCascade shapes and some return these wrappers. Therefore,
# we import all shape wrapper classes from the geometry module
from frackit.geometry import OCCShapeWrapper, OCCVertexWrapper, OCCEdgeWrapper, OCCWireWrapper
from frackit.geometry import OCCFaceWrapper, OCCShellWrapper, OCCSolidWrapper, OCCCompoundWrapper


###########################################################
# boolean operators for shape wrapper or geometry classes #
###########################################################

# helper function to get the wrapped shape of an object
def _getShapeWrapper(object):
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
    return _occutilities.cut(_getShapeWrapper(object), _getShapeWrapper(tool), eps)

# intersection operator
def intersect(object, tool, eps):
    """
    Boolean operation for intersections.

    Parameters:
    object: The shape whose intersection with the tool is to be computed.
    tool: The tool shape.
    eps: Tolerance value to be used.
    """
    return _occutilities.intersect(_getShapeWrapper(object), _getShapeWrapper(tool), eps)

# fragmentation operator
def fragment(shapeList, eps):
    """
    Boolean operation for fragmentation.

    Parameters:
    shapeList: List of shapes whose fragments are to be computed.
    eps: Tolerance value to be used.
    """
    return _occutilities.fragment([_getShapeWrapper(shape) for shape in shapeList], eps)

# fusion operator
def fuse(shapeList, eps):
    """
    Boolean operation for fusions.

    Parameters:
    shapeList: List of shapes whose fusion is to be computed.
    eps: Tolerance value to be used.
    """
    return _occutilities.fuse([_getShapeWrapper(shape) for shape in shapeList], eps)

############################
# Bounding box computation #
############################
def getBoundingBox(object):
    return _occutilities.getBoundingBox(_getShapeWrapper(object))

#################################
# Write function for geometries #
#################################
# function supports the file types supported by OpenCascade
def write(object, fileName):
    return _occutilities.write(_getShapeWrapper(object), fileName)
