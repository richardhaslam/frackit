from ._entitynetwork import *

class EntityNetworkConstraints(_entitynetwork._EntityNetworkConstraints):

    def evaluate(self, a, b):
        """
        Evaluates the constraints between the arguments a and b, both being
        either an instance of a geometry class or a list of geometry classes.

        Parameters:
        a: a geometry or list of geometries
        b: a geometry or list of geometries
        """
        isListA = isinstance(a, list)
        isListB = isinstance(b, list)

        def evalList(theList, geo):
            return all( super(EntityNetworkConstraints, self).evaluate(item, b) for item in theList )

        if isListA and not isListB:
            return evalList(a, b)
        elif isListB and not isListA:
            return evalList(b, a)
        elif isListA and isListB:
            return all( evalList(b, item) for item in a )
        # here we assume the two arguments to be geometries
        else:
            return super().evaluate(a, b)
